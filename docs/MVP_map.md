# Pipeline plan

## MVP — Arquitectura mínima viable

```
FASTQ / BAM input
    │
    ▼
STARsolo (CB+UB) ──► reference bundle (STAR idx + GTF + term. exons + PAS ref)
    │                       │
    ▼                       ▼
Barcode filter      PAS_REFERENCE_BUILD
    │               (annotation + PolyA_DB/PolyASite + merge@25bp)
    ▼
Grouped BAMs (cluster, cell_type)
    │
    ├──► Coverage (bedGraph + bigWig)
    │
    └──► PDUI comparison (APA_CALLING) ──► Internal Priming Filter
                                │
                                ▼
                         apa_events.tsv + pdui_usage_matrix.tsv
                                │
                                ▼
                           HTML Report
```

### Módulos core (siempre activos)

| Módulo | Descripción | Toggle |
|---|---|---|
| `STARSOLO_ALIGN_SC` | Alineación con extracción de CB+UB | — |
| `PREPARE_REFERENCE_BUNDLE` | Índice STAR + GTF + chrom.sizes | — |
| `PAS_REFERENCE_BUILD` | Referencia PAS desde anotación + atlas | `enable_pas_reference_build` |
| `BUILD_TERMINAL_EXON_CATALOG` | Exones terminales desde GTF | — |
| `BAM_BARCODE_FILTER` | Filtrado a barcodes celulares válidos | — |
| `GROUPED_BAM_GENERATION` | Split de BAM por cluster/cell_type | — |
| `COVERAGE_TRACKS` | bedGraph + bigWig strand-aware | — |
| `APA_FEATURE_EXTRACTION` | 8 features por sitio×grupo | — |
| `APA_CALLING` | PDUI + Mann-Whitney U entre grupos | — |
| `INTERNAL_PRIMING_FILTER` | Filtro blacklist BED (ligero) | — |
| `RENDER_APA_REPORT` | Reporte HTML interactivo | — |

---

## MVP extendido — Módulos opcionales

```
MVP +
    ├── SIERRA_QUANT       (cuantificación de isoformas PAS por grupo)
    ├── SCAPTURE           (modelo LSTM para filtrado de internal priming)
    ├── APA_MODEL_TRAIN / APA_MODEL_PREDICT   (sklearn sobre datos propios)
    ├── APARENT scores     (score de secuencia fusionado en feature vector)
    ├── Single-cell projection   (APA scores a nivel de célula individual)
    └── DEXSeq / GLM       (test diferencial por réplicas, R)
```

### Módulos opcionales y criterio de activación

| Módulo | Toggle | Cuándo activar |
|---|---|---|
| `SIERRA_QUANT` | `enable_sierra_quant` | Cuando haya wrapper real de Sierra R; actualmente scaffold stub |
| `PAS_SCORING` | `enable_pas_scoring` | Cuando haya lógica de scoring real; actualmente scaffold stub |
| `APA_MODEL_TRAIN/PREDICT` | siempre corre si hay eventos | Útil con labels supervisadas propias; output actual es descriptivo |
| `enable_single_cell_apa_projection` | param bool | Solo si los clusters tienen suficiente cobertura por célula |
| SCAPTURE | pendiente de integración | Reemplaza `INTERNAL_PRIMING_FILTER` cuando dep. de DL sea aceptable |
| APARENT | pendiente de integración | Agregar `aparent_score` como columna extra en `apa_features.tsv` |
| DEXSeq/GLM | pendiente de integración | Requiere R + réplicas biológicas reales |

---

## Secuencia de adopción recomendada

1. **Ahora**: validar MVP con pbmc1k — confirmar cobertura real, eventos APA, `pdui_usage_matrix.tsv` no vacío
2. **Siguiente**: promover `SIERRA_QUANT` de sidecar a default una vez implementado el wrapper R real
3. **Después**: integrar SCAPTURE como módulo opcional que puede reemplazar `INTERNAL_PRIMING_FILTER`
4. **Largo plazo**: APARENT score como feature adicional, DEXSeq/GLM para análisis diferencial con réplicas

