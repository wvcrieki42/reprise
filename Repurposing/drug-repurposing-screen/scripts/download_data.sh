#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Fetch the open bulk datasets needed for a full-scale run (all approved drugs).
# All sources are free; DrugBank is optional and requires a (free academic) licence.
# After download, run build_full_tables.py to materialise the canonical CSV/Parquet,
# then set mode: full in config.yaml and point paths at data/full/*.
# ---------------------------------------------------------------------------
set -euo pipefail
mkdir -p data/full && cd data/full

echo ">> ChEMBL (approved drugs, mechanisms, indications)  ~5 GB"
# Resolve the latest sqlite archive from the "latest" index.
CHEMBL_ARCHIVE="${CHEMBL_ARCHIVE:-}"
if [[ -z "${CHEMBL_ARCHIVE}" ]]; then
  CHEMBL_ARCHIVE="$(wget -qO- "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/" \
    | tr '"' '\n' \
    | grep -E '^chembl_[0-9]+_sqlite\.tar\.gz$' \
    | sort -V \
    | tail -n1 || true)"
fi
if [[ -z "${CHEMBL_ARCHIVE}" ]]; then
  echo "ERROR: Could not resolve latest ChEMBL sqlite archive from index." >&2
  exit 1
fi
wget -nc "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/${CHEMBL_ARCHIVE}"
CHEMBL_DIR="${CHEMBL_ARCHIVE%_sqlite.tar.gz}"
if [[ ! -d "${CHEMBL_DIR}" ]]; then
  tar -xzf "${CHEMBL_ARCHIVE}"
fi

echo ">> Open Targets parquet datasets (associations, targets, evidence, baselineExpression)"
# Browse https://platform.opentargets.org/downloads for the current FTP path/version.
OT_VER="${OT_VER:-24.06}"
OT_BASE_URL="https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/${OT_VER}/output/etl/parquet"

download_ot_parquet_dir() {
  local dataset="$1"
  local src="${OT_BASE_URL}/${dataset}/"
  local dst="${dataset}"
  local listing
  local root_links
  local subdir_links
  local sub_listing
  local sub_files
  local root_count=0
  local subdir_count=0
  local href
  local file
  local subdir

  mkdir -p "${dst}"
  listing="$(mktemp)"
  root_links="$(mktemp)"
  subdir_links="$(mktemp)"
  wget -qO "${listing}" "${src}"
  while IFS= read -r href; do
    if [[ "${href}" == \?* ]] || [[ "${href}" == /* ]]; then
      continue
    fi
    if [[ "${href}" == */ ]]; then
      printf '%s\n' "${href%/}" >> "${subdir_links}"
      subdir_count=$((subdir_count + 1))
      continue
    fi
    if [[ "${href}" =~ (_SUCCESS|\.parquet)$ ]]; then
      printf '%s\n' "${href}" >> "${root_links}"
      root_count=$((root_count + 1))
    fi
  done < <(
    grep -Eo 'href="[^\"]+"' "${listing}" \
      | cut -d'"' -f2 || true
  )
  rm -f "${listing}"

  if [[ "${root_count}" -eq 0 && "${subdir_count}" -eq 0 ]]; then
    echo "ERROR: No dataset links found for Open Targets dataset '${dataset}' at ${src}" >&2
    rm -f "${root_links}" "${subdir_links}"
    return 1
  fi

  echo "   - ${dataset}: root_files=${root_count}, subdirs=${subdir_count}"
  while IFS= read -r file; do
    [[ -z "${file}" ]] && continue
    wget -nc -P "${dst}" "${src}${file}"
  done < "${root_links}"

  while IFS= read -r subdir; do
    [[ -z "${subdir}" ]] && continue
    mkdir -p "${dst}/${subdir}"
    sub_listing="$(mktemp)"
    sub_files="$(mktemp)"
    wget -qO "${sub_listing}" "${src}${subdir}/"
    grep -Eo 'href="[^\"]+"' "${sub_listing}" \
      | cut -d'"' -f2 \
      | grep -Ev '^\?|^/|/$' \
      | grep -E '\.parquet$' > "${sub_files}" || true
    rm -f "${sub_listing}"
    while IFS= read -r file; do
      [[ -z "${file}" ]] && continue
      wget -nc -P "${dst}/${subdir}" "${src}${subdir}/${file}"
    done < "${sub_files}"
    rm -f "${sub_files}"
  done < "${subdir_links}"

  rm -f "${root_links}" "${subdir_links}"
}

download_ot_parquet_dir "associationByOverallDirect"
download_ot_parquet_dir "targets"
download_ot_parquet_dir "evidence"
download_ot_parquet_dir "baselineExpression"

echo ">> EFO ontology (disease parents, for novelty radius)"
wget -nc "http://www.ebi.ac.uk/efo/efo.obo"

echo ">> (Optional) STRING full human network for offline expansion"
wget -nc "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz" || true

echo ">> (Optional, licensed) DrugBank full database XML -> place full_database.xml here manually"

cat <<'EOF'

Done. Next:
  1) python scripts/build_full_tables.py     # writes data/full/*.csv (canonical schema)
  2) edit config.yaml: set 'mode: full' and repoint paths to data/full/*
  3) make full
EOF
