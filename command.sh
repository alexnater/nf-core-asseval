module load Java
export APPTAINER_CACHEDIR=$SCRATCH
export NXF_SINGULARITY_CACHEDIR=/data/users/anater/singularity_cache
export NXF_TEMP=$SCRATCH
nextflow run main.nf -profile unibe_ibu -params-file test/params.yaml
