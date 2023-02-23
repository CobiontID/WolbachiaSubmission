# Wolbachia Submission Pipeline

To run this Snakemake pipeline, you need to generate a config yaml file:

```
shortname: tolid
workingdirectory: directory in which the markerscan results are store
outdir: directory in which to store the final output files, e.g. /nfs/team135/wolbachia/
stsfile: /lustre/scratch123/tol/tolqc/track/tol_sts.tsv
scientific_name: full scientific name, e.g. Tromatobia_lineatoria
gfafile: location of the gfa file (needs to correspond with the assembly which markerscan used)
pwd_directory: directory in which to store temporary files generated by this pipeline
```
This pipeline can now be run with the following command:

```
snakemake --configfile $configfile --cores 10 --use-conda --conda-prefix $condadir -s $pipelinedir/Snakefile
```

Or with singularity:

```
singularity pull docker://emvcaest/wolbachia_submission:latest
singularity run wolbachia_submission.sif snakemake --cores $threads --use-conda --conda-prefix /opt/conda/ -s /WolbachiaSubmission/Snakefile --configfile $configfile
```