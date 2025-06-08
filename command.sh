module load Java
export APPTAINER_CACHEDIR=$SCRATCH
export NXF_SINGULARITY_CACHEDIR=/data/projects/p2025-0010_somatic_variant_calling_in_dog_glioblastoma/pipelines/singularity_cache
export NXF_TEMP=$SCRATCH
nextflow run main.nf -profile unibe_ibu -params-file test/params.yaml
