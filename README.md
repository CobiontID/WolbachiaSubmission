# Wolbachia Submission Pipeline

To run this Snakemake pipeline, you need to generate a config yaml file:

```
shortname: tolid
workingdirectory: directory in which the markerscan results are stored
outdir: directory in which to store the final output files, e.g. /lustre/scratch124/tol/projects/darwin/data/insects/$scientific_name/assembly/draft/
stsfile: /lustre/scratch123/tol/tolqc/track/tol_sts.tsv
scientific_name: full scientific name, e.g. Tromatobia_lineatoria
gfafile: location of the gfa file (needs to correspond with the assembly which markerscan used)
pwd_directory: directory in which to store temporary files generated by this pipeline
taxfile: location of the tax file from MarkerScan (*SSU.reduced.SILVA.tax)
email: email address from where to send emails to ENA taxid requests, if you leave this as an empty string emails will not be sent
```
You now can run this using Singularity:
```
singularity pull docker://{repo}/wolbachiasubmission:{tag}
singularity run wolbachiasubmission.sif snakemake --cores $threads --use-conda --conda-prefix /opt/conda/envs -s /WolbachiaSubmission/Snakefile --configfile $configfile
```
