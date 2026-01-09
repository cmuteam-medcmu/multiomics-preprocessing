#!/bin/bash

#SBATCH --job-name=dicom2nifti
#SBATCH --account=treephum
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=dicom2nifti%j.out
#SBATCH --error=dicom2nifti%j.err
#SBATCH --mail-type=END

# =================================================
# === 1. Environment Setup ===
# =================================================

micromamba activate dicom2nifti

# =================================================
# === 2. User Config (load from YAML) ===
# =================================================

CONFIG_FILE="/project/o250004_STSML01/Radiomics/dcm2niix/dcm2nii_config.yaml"

# ---- sanity checks ----
if [ ! -f "$CONFIG_FILE" ]; then
  echo "[ERROR] Missing config file: $CONFIG_FILE"
  exit 1
fi

if ! command -v yq &> /dev/null; then
  echo "[ERROR] yq not found. Install with: micromamba install -c conda-forge yq"
  exit 1
fi

sanitize_path() {
  echo "$1" | tr -d '\r' | tr -d '"' | xargs
}

# ---- params ----
ROOT=$(sanitize_path "$(yq -r '.params.ROOT' "$CONFIG_FILE")")
OUTDIR=$(sanitize_path "$(yq -r '.params.OUTDIR' "$CONFIG_FILE")")

MODALITY=$(sanitize_path "$(yq -r '.params.MODALITY' "$CONFIG_FILE")")
IMAGE_KEY=$(sanitize_path "$(yq -r '.params.IMAGE_KEY' "$CONFIG_FILE")")
MASK_KEY=$(sanitize_path "$(yq -r '.params.MASK_KEY' "$CONFIG_FILE")")
FIND_DEPTH=$(yq -r '.params.FIND_DEPTH' "$CONFIG_FILE")

mapfile -t CASE_PATTERNS < <(yq -r '.params.CASE_PATTERNS[]' "$CONFIG_FILE")

# optional: run tag prefix (default to "run" if missing)
RUN_TAG_PREFIX=$(sanitize_path "$(yq -r '.params.RUN_TAG_PREFIX // "run"' "$CONFIG_FILE")")

# ---- dcm2niix ----
DCM_GZ=$(sanitize_path "$(yq -r '.dcm2niix.gz' "$CONFIG_FILE")")
DCM_DEPTH=$(yq -r '.dcm2niix.depth' "$CONFIG_FILE")
DCM_QUICK=$(sanitize_path "$(yq -r '.dcm2niix.quick_search' "$CONFIG_FILE")")
DCM_IGNORE=$(sanitize_path "$(yq -r '.dcm2niix.ignore_derived' "$CONFIG_FILE")")
DCM_SINGLE=$(sanitize_path "$(yq -r '.dcm2niix.single' "$CONFIG_FILE")")
DCM_FILENAME=$(sanitize_path "$(yq -r '.dcm2niix.filename' "$CONFIG_FILE")")
DCM_WRITE=$(yq -r '.dcm2niix.write' "$CONFIG_FILE")
DCM_MERGE=$(yq -r '.dcm2niix.merge' "$CONFIG_FILE")
DCM_CROP=$(sanitize_path "$(yq -r '.dcm2niix.crop' "$CONFIG_FILE")")
DCM_SCALE=$(sanitize_path "$(yq -r '.dcm2niix.lossless_scale' "$CONFIG_FILE")")
DCM_PHILIPS=$(sanitize_path "$(yq -r '.dcm2niix.philips' "$CONFIG_FILE")")
DCM_ADJACENT=$(sanitize_path "$(yq -r '.dcm2niix.adjacent' "$CONFIG_FILE")")
DCM_RENAME=$(sanitize_path "$(yq -r '.dcm2niix.rename' "$CONFIG_FILE")")
DCM_EXPORT=$(sanitize_path "$(yq -r '.dcm2niix.export' "$CONFIG_FILE")")
DCM_VERBOSE=$(yq -r '.dcm2niix.verbose' "$CONFIG_FILE")

# ---- BIDS split ----
IMG_BIDS=$(sanitize_path "$(yq -r '.dcm2niix.bids_image' "$CONFIG_FILE")")
MSK_BIDS=$(sanitize_path "$(yq -r '.dcm2niix.bids_mask' "$CONFIG_FILE")")
ANON_BIDS=$(sanitize_path "$(yq -r '.dcm2niix.anonymize_bids' "$CONFIG_FILE")")

# ---- dcmrtstruct2nii ----
RT_COMMAND=$(sanitize_path "$(yq -r '.dcmrtstruct2nii.command' "$CONFIG_FILE")")
RT_MERGE=$(yq -r '.dcmrtstruct2nii.merge_rois' "$CONFIG_FILE")
RT_MERGE_LABEL=$(yq -r '.dcmrtstruct2nii.merged_label_value' "$CONFIG_FILE")
RT_STRICT=$(yq -r '.dcmrtstruct2nii.strict_geometry' "$CONFIG_FILE")
RT_FILL=$(yq -r '.dcmrtstruct2nii.fill_holes' "$CONFIG_FILE")
RT_ALLOW_EMPTY=$(yq -r '.dcmrtstruct2nii.allow_empty_masks' "$CONFIG_FILE")
RT_VERBOSE=$(yq -r '.dcmrtstruct2nii.verbose' "$CONFIG_FILE")

mapfile -t RT_INCLUDE < <(yq -r '.dcmrtstruct2nii.include_rois[]?' "$CONFIG_FILE")
mapfile -t RT_EXCLUDE < <(yq -r '.dcmrtstruct2nii.exclude_rois[]?' "$CONFIG_FILE")

# ---- field checks ----
for v in ROOT OUTDIR IMAGE_KEY MASK_KEY; do
  if [ -z "${!v}" ] || [ "${!v}" = "null" ]; then
    echo "[ERROR] Missing required config: $v"
    exit 1
  fi
done

if [ "${#CASE_PATTERNS[@]}" -eq 0 ]; then
  echo "[ERROR] params.CASE_PATTERNS is empty"
  exit 1
fi

# =================================================
# 3. Prepare output (per-run folder)
# =================================================

mkdir -p "$OUTDIR"

RUN_TAG="${RUN_TAG_PREFIX}_$(date +%Y%m%d_%H%M%S)"
RUN_DIR="$OUTDIR/$RUN_TAG"
mkdir -p "$RUN_DIR"
echo "[RUN_DIR] $RUN_DIR"

# ---- record the run config ----
cat > "$RUN_DIR/run_config.txt" <<EOF
CONFIG_FILE=$CONFIG_FILE

ROOT=$ROOT
OUTDIR=$OUTDIR
RUN_DIR=$RUN_DIR

CASE_PATTERNS=${CASE_PATTERNS[*]}
MODALITY=$MODALITY
IMAGE_KEY=$IMAGE_KEY
MASK_KEY=$MASK_KEY
FIND_DEPTH=$FIND_DEPTH

dcm2niix:
  gz=$DCM_GZ
  bids_image=$IMG_BIDS
  bids_mask=$MSK_BIDS
  anonymize_bids=$ANON_BIDS
  ignore_derived=$DCM_IGNORE
  filename=$DCM_FILENAME
  write=$DCM_WRITE
  merge=$DCM_MERGE
  crop=$DCM_CROP
  lossless_scale=$DCM_SCALE
  philips=$DCM_PHILIPS
  verbose=$DCM_VERBOSE

dcmrtstruct2nii:
  command=$RT_COMMAND
  strict_geometry=$RT_STRICT
  fill_holes=$RT_FILL
  allow_empty_masks=$RT_ALLOW_EMPTY
  verbose=$RT_VERBOSE
EOF

# =================================================
# 4. Auto detect CASE
# =================================================

CASE_DIRS=()

for pat in "${CASE_PATTERNS[@]}"; do
  while IFS= read -r d; do
    CASE_DIRS+=("$d")
  done < <(find "$ROOT" -mindepth 1 -maxdepth 2 -type d -name "$pat" 2>/dev/null)
done

IFS=$'\n' CASE_DIRS=($(printf "%s\n" "${CASE_DIRS[@]}" | sort -u))
unset IFS

if [ "${#CASE_DIRS[@]}" -eq 0 ]; then
  echo "[ERROR] No CASE directories detected under ROOT: $ROOT"
  exit 1
fi

echo "Detected CASE_DIRS:"
for c in "${CASE_DIRS[@]}"; do
  echo "$c"
done

# =================================================
# helpers (defined once, used for each case)
# =================================================

# ---- Find first directory under root (within maxdepth) whose basename matches regex (case-insensitive) ----
first_match_dir_under() {
  local root="$1"; local maxdepth="$2"; local re="$3"
  find "$root" -mindepth 1 -maxdepth "$maxdepth" -type d 2>/dev/null \
    | while IFS= read -r d; do
        local bn
        bn=$(basename "$d")
        echo "$bn" | grep -qiE "$re" && { echo "$d"; return 0; }
      done
  return 1
}

# ---- Pick first dir that matches include regex but NOT exclude regex (case-insensitive) ----
first_match_dir_under_excluding() {
  local root="$1"; local maxdepth="$2"; local include_re="$3"; local exclude_re="$4"

  find "$root" -mindepth 1 -maxdepth "$maxdepth" -type d 2>/dev/null \
    | while IFS= read -r d; do
        local bn
        bn=$(basename "$d")
        echo "$bn" | grep -qiE "$include_re" || continue
        if [ -n "$exclude_re" ]; then
          echo "$bn" | grep -qiE "$exclude_re" && continue
        fi
        echo "$d"
        return 0
      done
  return 1
}

# ---- List all directories under root (within maxdepth) whose basename matches regex (case-insensitive) ----
list_match_dirs_under() {
  local root="$1"; local maxdepth="$2"; local re="$3"
  find "$root" -mindepth 1 -maxdepth "$maxdepth" -type d 2>/dev/null \
    | while IFS= read -r d; do
        local bn
        bn=$(basename "$d")
        echo "$bn" | grep -qiE "$re" && echo "$d"
      done
}

# ---- Count regular files directly inside dir (depth 1) ----
count_files_shallow() {
  local dir="$1"
  [ -d "$dir" ] || { echo 0; return 0; }
  find "$dir" -mindepth 1 -maxdepth 1 -type f 2>/dev/null | wc -l
}

# ---- Pick first regular file directly inside dir (depth 1) ----
pick_first_file_shallow() {
  local dir="$1"
  [ -d "$dir" ] || return 1
  find "$dir" -mindepth 1 -maxdepth 1 -type f 2>/dev/null | head -n 1
}

# ---- Choose best first-level modality folder under CASE_DIR, excluding EX_RE. ----
pick_modality_dir_excluding() {
  # Prefer MRI/MR if present; otherwise pick non-excluded folder with most files (shallow).
  local case_dir="$1"
  local ex_re="$2"

  mapfile -t cands < <(find "$case_dir" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)
  [ "${#cands[@]}" -eq 0 ] && return 1

  # 1) Prefer MRI/MR (not excluded)
  for d in "${cands[@]}"; do
    local bn
    bn=$(basename "$d")
    echo "$bn" | grep -qiE "$ex_re" && continue
    echo "$bn" | grep -qiE '(^|[^a-z])mr(i)?([^a-z]|$)' && { echo "$d"; return 0; }
  done

  # 2) Otherwise pick non-excluded folder with most files (shallow)
  local best=""
  local best_n=-1
  for d in "${cands[@]}"; do
    local bn n
    bn=$(basename "$d")
    echo "$bn" | grep -qiE "$ex_re" && continue
    n=$(find "$d" -mindepth 1 -maxdepth 2 -type f 2>/dev/null | wc -l)
    if [ "$n" -gt "$best_n" ]; then
      best_n="$n"
      best="$d"
    fi
  done

  [ -n "$best" ] && { echo "$best"; return 0; }
  return 1
}

# =================================================
# MAIN LOOP: process each case
# =================================================
for CASE_DIR in "${CASE_DIRS[@]}"; do
  CASE=$(basename "$CASE_DIR")
  CASE_OUTDIR="$RUN_DIR/$CASE"
  mkdir -p "$CASE_OUTDIR"

  echo "========================================"
  echo "[CASE] $CASE_DIR"
  echo "[OUT ] $CASE_OUTDIR"
  echo "========================================"

  # =================================================
  # 5. Resolve WORK_DIR + IMGDIR + MASKDIRS
  # =================================================

  # 5.1 ---- resolve WORK_DIR (case-level or modality-level) ----
  WORK_DIR="$CASE_DIR"

  if [ -n "$MODALITY" ] && [ "$MODALITY" != "null" ]; then
    if [[ "$MODALITY" == -* ]]; then
      EX_RE_RAW="${MODALITY#-}"
      EX_RE="($EX_RE_RAW|pet|pet[_ -]?ct|ct)"
      MOD_DIR=$(pick_modality_dir_excluding "$CASE_DIR" "$EX_RE")
      [ -n "${MOD_DIR:-}" ] && WORK_DIR="$MOD_DIR"
    else
      MOD_DIR=$(first_match_dir_under "$CASE_DIR" 1 "$MODALITY")
      [ -n "${MOD_DIR:-}" ] && WORK_DIR="$MOD_DIR"
    fi
  fi

  echo "WORK_DIR : $WORK_DIR (MODALITY='$MODALITY')"

  # 5.2 ---- resolve IMGDIR ----
  IMG_EXCLUDE_RE="rtstruct|mask|seg|roi|$MASK_KEY"
  IMGDIR=$(first_match_dir_under_excluding "$WORK_DIR" "$FIND_DEPTH" "$IMAGE_KEY" "$IMG_EXCLUDE_RE")
  if [ -z "${IMGDIR:-}" ] || [ ! -d "$IMGDIR" ]; then
    echo "[SKIP] IMGDIR not found (IMAGE_KEY='$IMAGE_KEY') under: $WORK_DIR"
    continue
  fi
  echo "IMGDIR   : $IMGDIR (IMAGE_KEY='$IMAGE_KEY')"

  # 5.3 ---- resolve MASKDIRS (collect all) ----
  mapfile -t MASKDIRS < <(list_match_dirs_under "$WORK_DIR" "$FIND_DEPTH" "$MASK_KEY" | sort -u)
  echo "MASKDIRS (MASK_KEY='$MASK_KEY'):"
  if [ "${#MASKDIRS[@]}" -eq 0 ]; then
    echo "  NONE"
  else
    for md in "${MASKDIRS[@]}"; do
      echo "  $md"
    done
  fi

  # 5.4 ---- classify masks -> MASK_TASKS ----
  MASK_TASKS=()
  for md in "${MASKDIRS[@]}"; do
    nfiles=$(count_files_shallow "$md")
    if [ "$nfiles" -eq 1 ]; then
      rtfile=$(pick_first_file_shallow "$md")
      MASK_TASKS+=("RTSTRUCT|$md|$rtfile")
    else
      MASK_TASKS+=("SERIES|$md|")
    fi
  done

  echo "MASK_TASKS:"
  if [ "${#MASK_TASKS[@]}" -eq 0 ]; then
    echo "  NONE"
  else
    for t in "${MASK_TASKS[@]}"; do
      IFS='|' read -r typ mdir mfile <<< "$t"
      if [ "$typ" = "RTSTRUCT" ]; then
        echo "  [$typ] dir=$mdir file=$mfile  -> run 7.1"
      else
        echo "  [$typ] dir=$mdir              -> run 7.2"
      fi
    done
  fi

  # =================================================
  # 6. Convert image
  # =================================================
  echo "[RUN] dcm2niix image -> $CASE_OUTDIR"
  dcm2niix \
    -z "$DCM_GZ" \
    -b "$IMG_BIDS" \
    -i "$DCM_IGNORE" \
    -o "$CASE_OUTDIR" \
    "$IMGDIR"

  # =================================================
  # 7. Convert mask(s)
  # =================================================
  if [ "${#MASK_TASKS[@]}" -eq 0 ]; then
    echo "[INFO] No mask folders found"
  else
    for t in "${MASK_TASKS[@]}"; do
      IFS='|' read -r typ mdir mfile <<< "$t"
      bn=$(basename "$mdir")

      if [ "$typ" = "RTSTRUCT" ]; then
        # 7.1 ---- RTSTRUCT -> NIfTI ----
        echo "[RUN] 7.1 dcmrtstruct2nii : $mfile"
        dcmrtstruct2nii "$RT_COMMAND" \
          -r "$mfile" \
          -d "$IMGDIR" \
          -o "$CASE_OUTDIR" \
      else
        # 7.2 ---- Mask series -> NIfTI ----
        echo "[RUN] 7.2 dcm2niix mask series : $mdir"
        dcm2niix \
          -z "$DCM_GZ" \
          -b "$MSK_BIDS" \
          -i "$DCM_IGNORE" \
          -f "${CASE}_${bn}" \
          -o "$CASE_OUTDIR" \
          "$mdir"
      fi
    done
  fi

  echo "[DONE] $CASE"
  echo "----------------------------------------"
done