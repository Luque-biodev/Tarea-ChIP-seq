#!/usr/bin/env bash
set -euo pipefail

# ==========================================================
# download.sh
# Objetivo: descargar los picos desde GEO, dejarlos listos para
#           el script de R y ejecutar HOMER.
# Parámetros: se leen de params.txt (o de otro archivo pasado
#             como primer argumento: bash download.sh mis_params.txt).
# Requisitos: curl, gzip, sed y HOMER (configureHomer.pl y
#             findMotifsGenome.pl en el PATH).
# ==========================================================

# -----------------------------
# Cargar parámetros
# -----------------------------
PARAMS="${1:-params.txt}"
[[ -f "${PARAMS}" ]] || { echo "No se encuentra el archivo de parámetros: ${PARAMS}" >&2; exit 1; }
# shellcheck disable=SC1090
source "${PARAMS}"

# Nombres derivados de los parámetros
STEM="${GID}_${AGI}_${TF}_${CONDITION}_${PEAK_SUFFIX}"
OUT="${STEM}.bed.gz"

# GEO exige "_" y "." codificados en el parámetro file=
ENC="${OUT//_/%5F}"; ENC="${ENC//./%2E}"
URL="https://www.ncbi.nlm.nih.gov/geo/download/?acc=${GID}&format=file&file=${ENC}"

bed_gz="${DATA_DIR}/${OUT}"
bed_unz="${DATA_DIR}/${STEM}.bed"
# Mismo nombre que espera tarea_ChIP_seq.R
bed_nochr="${DATA_DIR}/${STEM}_reanotado.bed"
outdir="${RESULTS_DIR}/dnaMotifs"

# -----------------------------
# Comprobación de dependencias
# -----------------------------
for cmd in curl gzip sed configureHomer.pl findMotifsGenome.pl; do
    command -v "${cmd}" > /dev/null || { echo "Falta '${cmd}' en el PATH." >&2; exit 1; }
done

# -----------------------------
# Descarga (solo si no existe)
# -----------------------------
mkdir -p "${DATA_DIR}"
if [[ -f "${bed_gz}" ]]; then
    echo "El archivo ${bed_gz} ya existe; se omite la descarga."
else
    echo "Descargando ${OUT} desde GEO..."
    curl -L --fail -o "${bed_gz}.part" "${URL}"
    # Comprobar que es un gzip válido (GEO a veces devuelve HTML con código 200)
    gzip -t "${bed_gz}.part" || { rm -f "${bed_gz}.part"; echo "La descarga no es un .gz válido." >&2; exit 1; }
    mv "${bed_gz}.part" "${bed_gz}"
    echo "Descarga completa: ${bed_gz}"
fi

# -----------------------------
# Genoma de HOMER
# (se guarda la salida en una variable: con `set -o pipefail`,
#  `cmd | grep -q` puede fallar por SIGPIPE y reinstalar sin necesidad)
# -----------------------------
homer_list="$(configureHomer.pl -list || true)"
if grep -qE "^[[:space:]+-]*${GENOME}([[:space:]]|$)" <<< "${homer_list}"; then
    echo "Genoma '${GENOME}' ya instalado."
else
    echo "Genoma '${GENOME}' no encontrado. Instalando..."
    configureHomer.pl -install "${GENOME}"
fi

# -----------------------------
# Preparar el BED: descomprimir y quitar el prefijo 'chr'
# -----------------------------
gzip -dc "${bed_gz}" > "${bed_unz}"
sed 's/^chr//' "${bed_unz}" > "${bed_nochr}"

# -----------------------------
# HOMER
# -----------------------------
mkdir -p "${outdir}"
findMotifsGenome.pl "${bed_nochr}" "${GENOME}" "${outdir}" -size "${HOMER_SIZE}" -len "${HOMER_LEN}"

echo "HOMER completado. Resultados en ${outdir}"
if [[ -f "${outdir}/knownResults.txt" ]]; then
    head -n 20 "${outdir}/knownResults.txt"
else
    echo "No se encontró knownResults.txt en ${outdir}."
fi
echo
echo "Para el script de R, usa:  data_folder <- \"${DATA_DIR}\""
