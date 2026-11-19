# Flujo de análisis ChIP-seq: anotación de picos, genes diana y enriquecimiento

Entrega de los scripts, figuras y tablas del análisis de los picos de ChIP-seq de ABF3 (GSE80564, *Arabidopsis thaliana*). Este documento describe **cómo está organizado el flujo de trabajo y cómo reproducirlo**. La interpretación de los resultados está en el informe PDF.

## 1. Visión general del flujo

```
 params.txt ──► download.sh
                   │  [Bash] descarga desde GEO, descompresión y sed (quitar prefijo "chr")
                   ▼
        data/GSE80564/..._reanotado.bed
                   ├──────────────────────────────┐
                   ▼ [R] tarea_ChIP_seq.R         ▼ [Bash] download.sh (continuación)
             ChIPseeker: annotatePeak          HOMER: findMotifsGenome.pl
                   │                           (known + de novo)
                   ├─ QC: pie, bar, distancia a TSS
                   ├─ anotación manual (bucle propio)
                   ▼
             genes diana (picos en promotor)
                   ▼
             enrichGO (BP)  ·  enrichKEGG  ·  pathview
```

**Alcance:** el flujo parte de los **picos ya llamados** (archivo `optimal_narrowPeak` publicado en GEO). Los pasos previos (FastQC, Bowtie2, Samtools, MACS2) los hicieron los autores del estudio original y **no se reejecutan aquí**.

## 2. Contenido de la entrega

| Archivo | Tipo | Papel en el flujo |
|---|---|---|
| `params.txt` | Configuración | Parámetros de la descarga y de HOMER (serie GEO, factor, genoma, `-size`, `-len`, rutas). |
| `download.sh` | Script Bash | Lee `params.txt`, descarga los picos de GEO, los prepara y ejecuta HOMER. |
| `tarea_ChIP_seq.R` | Script R | Anotación, genes diana, GO y KEGG. |
| `tarea_ChIP_seq3.Rmd` | Informe reproducible | El mismo código de R, con texto y figuras. |
| `tarea_ChIP_seq_ABF3_Alicia_y_Pablo.pdf` | Informe | Resultado de renderizar el `.Rmd`; incluye figuras que no salen de los scripts (IGV, HOMER). |
| `target_genes_from_ABF3_peaks.txt` | Salida | Genes diana, un ID TAIR por línea (9 400). |
| `KEGG_ABF3_enrichment.tsv` | Salida | Tabla KEGG (en realidad CSV, por usar `write.csv`). |
| `ABF3_peak_annotation_pie.png`, `ABF3_peak_dist_to_TSS.png` | Salida (QC) | Anotación y distancia al TSS (ChIPseeker). |
| `freq_peak_annotation_pie.png` | Salida (QC) | Anotación con la función manual. |
| `GO_ABF3_barplot.png`, `_dotplot.png`, `_emapplot.png`, `_cnetplot.png` | Salida | Gráficos de enriquecimiento GO. |

**Estructura de carpetas tras ejecutar todo** (se crea automáticamente):

```
.
├── params.txt
├── download.sh
├── tarea_ChIP_seq.R
├── data/GSE80564/
│   ├── GSE80564_AT4G34000_ABF3_ABA_optimal_narrowPeak_p16.bed.gz
│   ├── GSE80564_AT4G34000_ABF3_ABA_optimal_narrowPeak_p16.bed
│   └── GSE80564_AT4G34000_ABF3_ABA_optimal_narrowPeak_p16_reanotado.bed
└── resultsABF3/
    ├── dnaMotifs/        (salida de HOMER)
    └── ...               (figuras y tablas de R)
```

## 3. Parámetros, entradas y requisitos

### Archivo `params.txt`

Contiene los parámetros del paso de Bash, con formato `CLAVE="valor"`. `download.sh` lo carga con `source`.

| Parámetro | Valor por defecto | Significado |
|---|---|---|
| `GID` | `GSE80564` | Serie GEO. |
| `AGI` | `AT4G34000` | Identificador AGI del factor de transcripción. |
| `TF` | `ABF3` | Nombre del factor (también da nombre a la carpeta de resultados). |
| `CONDITION` | `ABA` | Tratamiento. |
| `PEAK_SUFFIX` | `optimal_narrowPeak_p16` | Tipo de archivo de picos publicado en GEO. |
| `GENOME` | `tair10` | Genoma de HOMER. |
| `HOMER_SIZE` | `100` | `-size` de HOMER: ventana total de 100 pb centrada en el pico (±50 pb). |
| `HOMER_LEN` | `8` | `-len` de HOMER: longitud del motivo. |
| `DATA_DIR` | `data/${GID}` | Carpeta de los datos descargados. |
| `RESULTS_DIR` | `results${TF}` | Carpeta de salidas, la misma que usa el script de R. |

Para usar otro archivo de parámetros: `bash download.sh mis_params.txt`.

### Parámetros del script de R

Estos están fijados dentro de `tarea_ChIP_seq.R` y no se leen de `params.txt`:

| Paso | Parámetro | Valor |
|---|---|---|
| Anotación de picos | `tssRegion` | `c(-1000, 1000)` |
| Definición de genes diana | criterio | `annotation == "Promoter"` |
| GO | `ont`, `keyType`, `pvalueCutoff`, `qvalueCutoff` | `BP`, `TAIR`, 0,05, 0,05 |
| KEGG | `organism`, `pAdjustMethod`, `pvalueCutoff` | `ath`, `BH`, 0,05 |

### Entrada

`GSE80564_AT4G34000_ABF3_ABA_optimal_narrowPeak_p16.bed.gz` de GEO. **`download.sh` la descarga automáticamente**, si el archivo no existe ya en `data/GSE80564/`.

### Software

Versiones tomadas del log del informe.

- Herramientas de línea de comandos: `curl`, `gzip`, `sed`.
- **HOMER** instalado, con `configureHomer.pl` y `findMotifsGenome.pl` en el `PATH`. El script no instala HOMER, pero sí descarga el genoma `tair10` si falta.
- R con Bioconductor: `ChIPseeker` 1.34.1, `clusterProfiler` 4.14.4, `enrichplot` 1.26.5, `TxDb.Athaliana.BioMart.plantsmart28`, `org.At.tair.db`, `pathview`.
- Conexión a internet: GEO (descarga), HOMER (genoma) y KEGG (`enrichKEGG`, `pathview`).

Para dejar constancia exacta del entorno, ejecuta al final del script de R `sessionInfo()` y guarda la salida.

## 4. Cómo reproducirlo, paso a paso

Ejecuta todo **desde la carpeta raíz del proyecto**, porque los scripts usan rutas relativas.

1. **Revisar `params.txt`.** Los valores por defecto reproducen el análisis de ABF3. Solo hay que editarlo para analizar otro factor.
2. **Descargar y preparar los datos, y ejecutar HOMER (Bash).**
   ```bash
   bash download.sh
   ```
   El script comprueba que estén `curl`, `gzip`, `sed` y HOMER en el `PATH` y después:
   - descarga el `.bed.gz` y comprueba que sea un gzip válido,
   - lo descomprime y quita el prefijo `chr` con `sed`, generando `..._reanotado.bed`,
   - instala el genoma `tair10` si no está,
   - ejecuta `findMotifsGenome.pl` con `HOMER_SIZE` y `HOMER_LEN`, con salida en `resultsABF3/dnaMotifs/`.
3. **Apuntar el script de R a los datos.** En `tarea_ChIP_seq.R`, cambia
   ```r
   data_folder <- file.path(".")
   ```
   por
   ```r
   data_folder <- file.path("data", "GSE80564")
   ```
4. **Ejecutar el análisis en R:** `tarea_ChIP_seq.R` (o *Knit* del `.Rmd`, con el mismo cambio de ruta). Realiza, por orden:
   - anotación de los picos (`readPeakFile` + `annotatePeak`),
   - control de calidad (`plotAnnoPie`, `plotAnnoBar`, `plotDistToTSS`) y anotación manual propia, con tabla en `ABF3.table.peak.annotation.tsv`,
   - definición del reguloma (picos en promotor) y guardado de `target_genes_from_ABF3.peaks.txt`,
   - enriquecimiento `enrichGO` y `enrichKEGG`, con sus tablas y gráficos; `pathview` dibuja una ruta KEGG (la 8.ª fila de la tabla).
5. **Revisar los motivos.** Los resultados de HOMER están en `resultsABF3/dnaMotifs/` (`knownResults.html`, `homerResults.html`). Las tablas del informe (figuras 12 y 13) son capturas de esos archivos.

El PDF se obtiene con *Knit* en RStudio sobre `tarea_ChIP_seq3.Rmd`. Además de los datos, el informe necesita las imágenes de IGV y HOMER, que se insertan a mano.

## 5. Limitaciones de reproducibilidad

- **`params.txt` solo configura el paso de Bash.** El script de R tiene sus parámetros y nombres (`ABF3`, rutas) fijados en el código. Si cambias `AGI` o `TF`, hay que adaptar también ese script.
- **La descarga depende de GEO.** Si el enlace cambia o GEO no responde, `download.sh` falla. Si ya tienes el `.bed.gz` en `data/GSE80564/`, se omite la descarga.
- **HOMER no se instala desde el script.** Hay que tenerlo configurado antes.
- **`enrichKEGG` y `pathview` leen KEGG en línea.** Si la base cambia, las tablas pueden variar ligeramente con respecto a las de la entrega.
- **Los números de la entrega (p. ej. 11 273 picos, 9 400 genes) corresponden a las versiones de arriba.** Con otras versiones de las anotaciones pueden diferir.
- **Hay un cambio manual de ruta** (`data_folder`) entre el script de Bash y el de R, y las capturas de IGV y de los motivos del informe se insertan a mano.
- **`KEGG_ABF3_enrichment.tsv` es CSV** aunque tenga extensión `.tsv`. Léelo con `read.csv()`.
- **El script de R no fija semilla ni registra el entorno.** Los pasos usados son deterministas, pero conviene añadir `sessionInfo()`.
- **El bucle manual de anotación recorre todos los picos con `subset()`** y tarda en ejecutarse con 11 273 picos.

## Referencias

- Song et al. (2016). *Science*, 354(6312), aag1550.
