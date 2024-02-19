"""
Pipeline to detect assembled Wolbachia by MarkerScan
----------------------------------------------------
Requirements 
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)
Basic usage:
  snakemake -p --use-conda --conda-prefix condadir --configfile config.yaml
"""
scriptdir = workflow.basedir+"/Scripts"
SSUHMMfile = workflow.basedir+"/SSU_Prok_Euk_Microsporidia.hmm"
SSUalnfile = workflow.basedir+"/SSU.aln.fa"
ncbi_taxdir = workflow.basedir+"/taxonomy"

pwd = config["workingdirectory"]
pwd_dir = config["pwd_directory"]
tolid = config["shortname"]
host = config["scientific_name"]
outdir = config["outdir"]
stsfile = config["stsfile"]
gfa = config["gfafile"]
email = config["email"]
taxreq = config["taxreq"]

rule all:
	input:
		expand('{pwd_dir}/assembly_done.txt', pwd_dir=config["pwd_directory"]),
		expand('{pwd_dir}/buscoAssembly/completeness_per_contig.txt', pwd_dir=config["pwd_directory"]),
		expand('{pwd_dir}/bins_done.txt', pwd_dir=config["pwd_directory"]),
		expand('{pwd_dir}/bins_SSU_done.txt', pwd_dir=config["pwd_directory"]),
		expand('{pwd_dir}/circular_done.txt', pwd_dir=config["pwd_directory"]),
		expand('{pwd_dir}/linear_done.txt', pwd_dir=config["pwd_directory"]),
		expand('{pwd_dir}/mags_done.txt', pwd_dir=config["pwd_directory"]),
		expand('{pwd_dir}/Wolbachia.done.txt', pwd_dir=config["pwd_directory"])


rule CheckPresenceWolbachia:
	"""
	check if this sample contains Wolbachia
	"""
	input:
		taxfile = expand("{tax}", tax=config["taxfile"])
	output:
		wolpresence = "{pwd_directory}/Wol_presence.txt"
	shell:
			"""
			if grep -q Wolbachia {input.taxfile} && [ -d {pwd} ]; then
				touch {output.wolpresence}
			fi
			"""

rule RunHifiasm:
	"""
	check if hifiasm was already run, if not re-try
	"""  
	input:
		wolpresence = "{pwd_directory}/Wol_presence.txt",
		fasta_reads = expand("{pwd}/kraken.fa", pwd=config["workingdirectory"])
	params:
		assemblyprefix = "{pwd_directory}/hifiasm/hifiasm",
		dirname = directory(expand("{pwd}/hifiasm",pwd=config["workingdirectory"])),
		gfa2 = expand("{pwd}/hifiasm/hifiasm.p_ctg.noseq.gfa",pwd=config["workingdirectory"])
	output:
		completed = "{pwd_directory}/assembly_done.txt",
		dirname2 = directory("{pwd_directory}/hifiasm"),
		gfa = "{pwd_directory}/hifiasm/hifiasm.p_ctg.gfa",
		fasta = "{pwd_directory}/hifiasm/hifiasm.p_ctg.fasta",
		fai = "{pwd_directory}/hifiasm/hifiasm.p_ctg.fasta.fai"
	threads: 10
	conda: "envs/hifiasm_seqtk.yaml"
	shell:
            """
			if [ ! -d {output.dirname2} ]; then
				mkdir {output.dirname2}
				cp {input.fasta_reads} {pwd_dir}/hifiasm/reads.fa
			fi
            if [ -s {input.fasta_reads} ] && [ ! -s {params.gfa2} ] ; then
				echo {input.fasta_reads}
				linecount=$(grep -c '>' < {input.fasta_reads})
				if [ $linecount -ge 50000 ]; then
					seqtk sample {input.fasta_reads} 50000 > {pwd_dir}/hifiasm/reads.fa
				fi
				hifiasm -o {params.assemblyprefix} -t {threads} {pwd_dir}/hifiasm/reads.fa -D 10 -l 1 -s 0.999 || true
				awk '/^S/{{print ">"$2"\\n"$3}}' {output.gfa} | fold > {output.fasta} || true
				faidx {output.fasta}
			else
				cp {params.dirname}/hifiasm.p_ctg.gfa {output.gfa} 
				cp {params.dirname}/hifiasm.p_ctg.fasta {output.fasta} 
				cp {params.dirname}/hifiasm.p_ctg.fasta.fai {output.fai}
			fi
			touch {output.completed}
			"""

rule RunBuscoAssembly:
    """
	Detect number of BUSCO genes per contig
	"""
	input:
		hifiasm_done = "{pwd_directory}/assembly_done.txt",
		circgenome = "{pwd_directory}/hifiasm/hifiasm.p_ctg.fasta",
	params:
		buscodir = directory(expand("{pwd}/buscoAssembly",pwd=config["workingdirectory"])),
		buscodir2 = directory("{pwd_directory}/buscoAssembly"),
		buscodbs = "{pwd_directory}/info_dbs_assembly.txt",
		buscoini = "{pwd_directory}/config_busco_assembly.ini"
	output:
		completed = "{pwd_directory}/buscoAssembly.done.txt",
		summary = "{pwd_directory}/buscoAssembly/completeness_per_contig.txt"
	conda: "envs/busco.yaml"
	threads: 10
	shell:
            """
			if [ ! -d {params.buscodir2} ]; then
				mkdir {params.buscodir2}
			fi

			if [ -f {params.buscodir}/busco/run*/full_table.tsv ] || [ -f {params.buscodir}/full_table.tsv ]; then
				if [ -f {params.buscodir}/full_table.tsv ]; then
					cp {params.buscodir}/full_table.tsv {params.buscodir2}
				else
					cp {params.buscodir}/busco/run*/full_table.tsv {params.buscodir2}
				fi
			else
				if [ -s {input.circgenome} ] ; then
					busco --list-datasets > {params.buscodbs}
                	python {scriptdir}/BuscoConfig.py -f {input.circgenome} -d {params.buscodir2} -dl busco_data/ -c {threads} -db {params.buscodbs} -o {params.buscoini}
                	busco --config {params.buscoini} -f || true
					mv {params.buscodir2}/busco/run*/full_table.tsv {params.buscodir2}
				else
					touch {params.buscodir2}/full_table.tsv
				fi
			fi
            touch {output.completed}
			python {scriptdir}/ParseBuscoTableMapping.py -d {params.buscodir2} -i {input.circgenome} -o {output.summary} 
            """

rule RunBusco:
    """
	Detect number of BUSCO genes per contig
	"""
	input:
		hifiasm_done = "{pwd_directory}/assembly_done.txt",
	params:
		circgenome =  expand("{pwd}/Anaplasmataceae.finalassembly.fa",pwd=config["workingdirectory"]),
		buscodir = directory(expand("{pwd}/busco",pwd=config["workingdirectory"])),
		buscodir2 = directory("{pwd_directory}/busco"),
		buscodbs = "{pwd_directory}/info_dbs_assembly.txt",
		buscoini = "{pwd_directory}/config_busco.ini",
	output:
		completed = "{pwd_directory}/busco.done.txt",
		summary = "{pwd_directory}/busco/completeness_per_contig.txt"
	conda: "envs/busco.yaml"
	threads: 10
	shell:
            """
			if [ ! -d {params.buscodir2} ]; then
				mkdir {params.buscodir2}
			fi

			if [ -f {params.buscodir}/busco/run*/full_table.tsv ] || [ -f {params.buscodir}/full_table.tsv ]; then
				if [ -f {params.buscodir}/full_table.tsv ]; then
					cp {params.buscodir}/full_table.tsv {params.buscodir2}
				else
					cp {params.buscodir}/busco/run*/full_table.tsv {params.buscodir2}
				fi
			else
				if [ -s {params.circgenome} ] ; then
					busco --list-datasets > {params.buscodbs}
                	python {scriptdir}/BuscoConfig.py -f {params.circgenome} -d {params.buscodir2} -dl busco_data/ -c {threads} -db {params.buscodbs} -o {params.buscoini}
                	busco --config {params.buscoini} -f || true
					mv {params.buscodir2}/busco/run*/full_table.tsv {params.buscodir2}
				else
					touch {params.buscodir2}/full_table.tsv
				fi
			fi
            touch {output.completed}
			python {scriptdir}/ParseBuscoTableMapping.py -d {params.buscodir2} -i {params.circgenome} -o {output.summary} 
            """

rule WolbachiaQuality:
	"""
	check quality of assemblies, create metadata file and chrom list file (if necessary). if quality not sufficient, create mag files
	"""
	input:
		hifiasm_done = "{pwd_directory}/assembly_done.txt",
		busco_asm_done = "{pwd_directory}/buscoAssembly.done.txt",
		busco_done = "{pwd_directory}/busco.done.txt",
		summary = "{pwd_directory}/buscoAssembly/completeness_per_contig.txt",
		orig_dir = directory(expand("{pwd}/",pwd=config["workingdirectory"]))
	output:
		completed = "{pwd_directory}/bins_done.txt",
		binlist = "{pwd_directory}/bin_list.txt"
	shell:
		"""
		python {scriptdir}/WolbachiaQuality.py -d {input.summary} -t {tolid} -l {output.binlist} -g {gfa} -d2 {input.orig_dir} -o {outdir}
		touch {output.completed}
		"""

checkpoint Bins:
	"""
	Move all fasta bins to directory supergroup
	"""
	input:
		bins_quality = "{pwd_directory}/bins_done.txt",
		binlist = "{pwd_directory}/bin_list.txt"
	output:
		superdir = directory("{pwd_directory}/supergroup/")
	shell:
		"""
		if [ ! -d {output.superdir} ]; then
			mkdir {output.superdir}
        fi
 		while read p
		do
			shortname=`echo $p | cut -d, -f1`
			echo {output.superdir}/$shortname.fa
			if [ ! -f {output.superdir}/bin.$shortname.fa ]; then
				mv {pwd_dir}/$shortname.fa {output.superdir}/bin.$shortname.fa
			fi
        done < {input.binlist}
		"""

rule HMMscan_SSU:
	"""
	Run HMMscan with prokaryotic+viral HMM (RF00177+RF01959)
	"""
	input:
		fasta = "{pwd_directory}/supergroup/bin.{bin}.fa"
	output:
		domfile = "{pwd_directory}/supergroup/{bin}.ProkSSU.domout",
		log = "{pwd_directory}/supergroup/{bin}.HMMscan.log"
	threads: 10
	conda: "envs/hmmer.yaml"
	shell:
		"""
		if [ -s {input.fasta} ]; then
			nhmmscan --cpu {threads} --noali --tblout {output.domfile} -o {output.log} {SSUHMMfile} {input.fasta}
		else
			touch {output.domfile} {output.log}
		fi
		"""

rule Fetch16SLoci:
	"""
	Get fasta sequences for detected reads with prokaryotic 16S signature and extract 16S locus
	"""
	input:
		dom = "{pwd_directory}/supergroup/{bin}.ProkSSU.domout", 
		fasta = "{pwd_directory}/supergroup/bin.{bin}.fa"
	output:
		readsinfo = "{pwd_directory}/supergroup/{bin}.ProkSSU.readsinfo",
		fasta16SLoci = "{pwd_directory}/supergroup/{bin}.ProkSSU.fa",
	conda: "envs/cdhit.yaml"
	shell:
		"""
		if [ -s {input.fasta} ]; then
			python {scriptdir}/GetReadsSSU_nhmmscan.py -i {input.dom} | grep -v 'RF02542.afa' > {output.readsinfo}
			python {scriptdir}/FetchSSUReads.py -i {output.readsinfo} -f {input.fasta} -o {output.fasta16SLoci}
		else
			touch {output.readsinfo} {output.fasta16SLoci}
		fi
		"""

rule Build_SSU_Tree:
	"""
	Add 16S to Wolbachia SSU alignment and build phylogenetic tree
	"""
	input:
		fasta16SLoci = "{pwd_directory}/supergroup/{bin}.ProkSSU.fa"
	output:
		aln16SLoci = "{pwd_directory}/supergroup/{bin}.ProkSSU.aln.fa",
		treefile = "{pwd_directory}/supergroup/{bin}.ProkSSU.treefile"
	conda: "envs/treebuild.yaml"
	threads: 10
	shell:
		"""
		if [ -s {input.fasta16SLoci} ]; then
			mafft --add {input.fasta16SLoci} --reorder {SSUalnfile} > {output.aln16SLoci}
			iqtree -s {output.aln16SLoci} -T {threads} -B 5000 -o Anaplasma_marginale_AF414871.1,Ehrlichia_canis_KJ513196.1
			mv {output.aln16SLoci}.treefile {output.treefile}
			rm {output.aln16SLoci}.*
		else
			touch {output.aln16SLoci} {output.treefile}
		fi
		"""
rule DownloadNCBITaxonomy:
        """
        Download current version of NCBI taxonomy
        """
        output:
                donefile = temporary("{pwd_directory}/taxdownload.done.txt")
        shell:
                """
                if [ ! -d {ncbi_taxdir} ]; then
                        mkdir {ncbi_taxdir}
                fi
                if [ -s {ncbi_taxdir}/names.dmp ]; then
                        before=$(date -d 'today - 180 days' +%s)
                        timestamp=$(stat -c %y {ncbi_taxdir}/names.dmp | cut -f1 -d ' ')
                        timestampdate=$(date -d $timestamp +%s)
                fi
        if [ ! -s {ncbi_taxdir}/names.dmp ] || [ $before -ge $timestampdate ]; then
                        curl -R https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz.md5 --output {ncbi_taxdir}/taxdump.tar.gz.md5
                        curl -R https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz --output taxdump.tar.gz
                        tar -C {ncbi_taxdir} -xzf taxdump.tar.gz names.dmp nodes.dmp
			rm taxdump.tar.gz
                fi
                touch {output.donefile}
                """

rule Supergroup:
	"""
	Link novel Wolbachia SSU to a supergroup
	"""
	input:
		treefile = "{pwd_directory}/supergroup/{bin}.ProkSSU.treefile",
		fasta16SLoci = "{pwd_directory}/supergroup/{bin}.ProkSSU.fa",
		taxfile = expand("{tax}", tax=config["taxfile"]),
		taxdone = "{pwd_directory}/taxdownload.done.txt"
	params:
		binname = "{bin}"
	output:
		supergroup_name =  "{pwd_directory}/supergroup/{bin}.Wol.SSU.supergroup.txt"
	conda: "envs/ete3.yaml"
	shell:	
		"""
		if [ -s {input.fasta16SLoci} ]; then
			python {scriptdir}/WolbachiaSupergroup.py -t {input.treefile} -c {input.fasta16SLoci} -b {params.binname} -s {host} -o {output.supergroup_name} -sts {stsfile} -i {tolid} -na {ncbi_taxdir}/names.dmp -no {ncbi_taxdir}/nodes.dmp -ta {input.taxfile}
			speciesname=` cat {output.supergroup_name} | cut -f1`
			echo "species: "$speciesname >>  {outdir}/{params.binname}*/{params.binname}*.yaml
			#cp {output.supergroup_name} {outdir}/{params.binname}.speciesname.txt
			#write code here to add novel name to file
			cat {output.supergroup_name} >> {pwd_dir}/taxid.nr.txt
		else
			touch {output.supergroup_name}
		fi
		"""

def aggregate_bins(wildcards):
	checkpoint_output=checkpoints.Bins.get(**wildcards).output[0]
	return expand(os.path.join(checkpoint_output, "{bin}.Wol.SSU.supergroup.txt"), bin=glob_wildcards(os.path.join(checkpoint_output, 'bin.{bin}.fa')).bin)

rule concatenate_bins:
	input:
		aggregate_bins
	output:
		"{pwd_directory}/bins_SSU_done.txt"
	shell:
		"""
		#write code here to submit email to ENA taxid request
		cut -f1-5 {pwd_dir}/taxid.nr.txt | sort | uniq > {pwd_dir}/taxid_request.txt

		#add script double-check it doesn't exist
		python {scriptdir}/taxid_requests.py -i {pwd_dir}/taxid_request.txt -o {pwd_dir}/taxid_request.novel.txt
		sh {scriptdir}/email_request_taxids.sh {pwd_dir}/taxid_request.novel.txt {email}

		rm {pwd_dir}/taxid.nr.txt
		touch {output}
		"""

checkpoint Bin_Circular:
	"""
	Move all fasta bins to directory supergroup
	"""
	input:
		binsdone = "{pwd_directory}/bins_SSU_done.txt",
		bins_quality = "{pwd_directory}/bins_done.txt",
		binlist = "{pwd_directory}/bin_list.txt"
	output:
		superdir = directory("{pwd_directory}/circular/")
	shell:
		"""
		if [ ! -d {output.superdir} ]; then
			mkdir {output.superdir}
        fi
 		while read p
		do
			shortname=`echo $p | cut -d, -f1`
			echo {pwd_dir}/$shortname.list
			if [ -f {pwd_dir}/$shortname.list ] && grep "Circular" {pwd_dir}/$shortname.list ; then
				echo {output.superdir}/bin.$shortname.fa
				touch {output.superdir}/bin.$shortname.txt
			fi
        done < {input.binlist}
		"""

rule AnnotateRotate:
	"""
	Annotate circular genomes to rotate at HemE
	"""
	input:
		fasta = "{pwd_directory}/supergroup/bin.{circ}.fa",
		binname= "{pwd_directory}/circular/bin.{circ}.txt"
	output:
		circ_fa = "{pwd_directory}/circular/{circ}.fa"
	params:
		circname = "{circ}"
	threads: 10
	conda: "envs/prokka.yaml"
	shell:
		"""
		prokka --cpus {threads} --outdir {pwd_dir}/circular/{params.circname}.prokka --prefix {params.circname} {input.fasta}
		python {scriptdir}/RotateRevComp.py -gff {pwd_dir}/circular/{params.circname}.prokka/{params.circname}.gff -fa {input.fasta} -o {output.circ_fa}
		today="$(date +'%Y%m%d')"
		rm {outdir}"/"{params.circname}"."$today"/"{params.circname}".fa.gz"
		cp {output.circ_fa} {outdir}"/"{params.circname}"."$today
		gzip {outdir}"/"{params.circname}"."$today"/"{params.circname}".fa"

		sh cobiont_curation_request.sh {outdir}"/"{params.circname}"."$today"/"{params.circname}".yaml" {email}
		"""

def aggregate_circs(wildcards):
	checkpoint_output=checkpoints.Bin_Circular.get(**wildcards).output[0]
	return expand(os.path.join(checkpoint_output, "{circ}.fa"), circ=glob_wildcards(os.path.join(checkpoint_output, 'bin.{circ}.txt')).circ)

rule concatenate_circular:
	input:
		aggregate_circs
	output:
		"{pwd_directory}/circular_done.txt"
	shell:
		"""
		touch {output}
		"""

checkpoint Bin_Linear:
	"""
	Move all fasta bins to directory supergroup
	"""
	input:
		binsdone = "{pwd_directory}/bins_SSU_done.txt",
		bins_quality = "{pwd_directory}/bins_done.txt",
		binlist = "{pwd_directory}/bin_list.txt"
	output:
		superdir = directory("{pwd_directory}/linear/")
	shell:
		"""
		if [ ! -d {output.superdir} ]; then
			mkdir {output.superdir}
        fi
 		while read p
		do
			shortname=`echo $p | cut -d, -f1`
			if [ -f {pwd_dir}/$shortname.list ] && grep "Linear" {pwd_dir}/$shortname.list ; then
				touch {output.superdir}/bin.$shortname.txt
				sh cobiont_curation_request.sh {outdir}"/"{params.circname}"."$today"/"{params.circname}".yaml" {email}
			elif [ ! -f {pwd_dir}/$shortname.list ] && [ ! -f {pwd_dir}/$shortname.manifest.txt ] && [ -f {pwd_dir}/$shortname.metadata.txt ] ; then
				touch {output.superdir}/bin.$shortname.txt
				sh cobiont_curation_request.sh {outdir}"/"{params.circname}"."$today"/"{params.circname}".yaml" {email}
			fi
        done < {input.binlist}
		"""

def aggregate_lins(wildcards):
	checkpoint_output=checkpoints.Bin_Linear.get(**wildcards).output[0]
	return expand(os.path.join(checkpoint_output, "bin.{circ}.txt"), circ=glob_wildcards(os.path.join(checkpoint_output, 'bin.{circ}.txt')).circ)

rule concatenate_linear:
	input:
		aggregate_lins
	output:
		"{pwd_directory}/linear_done.txt"
	shell:
		"""
		touch {output}
		"""

checkpoint Bin_Mags:
	"""
	Move all fasta bins to directory supergroup
	"""
	input:
		binsdone = "{pwd_directory}/bins_SSU_done.txt",
		bins_quality = "{pwd_directory}/bins_done.txt",
		binlist = "{pwd_directory}/bin_list.txt",
	output:
		superdir = directory("{pwd_directory}/mag/")
	shell:
		"""
		if [ ! -d {output.superdir} ]; then
			mkdir {output.superdir}
        fi
 		while read p
		do
			shortname=`echo $p | cut -d, -f1`
			if [ -f {pwd_dir}/$shortname.manifest.txt ] ; then
				mv {pwd_dir}/$shortname.manifest.txt {output.superdir}
				#cp {output.superdir}/$shortname.manifest.txt {outdir}
				#cp {pwd_dir}/supergroup/bin.$shortname.fa {outdir}/$shortname.mag.fa 
			fi
        done < {input.binlist}
		"""

def aggregate_mags(wildcards):
	checkpoint_output=checkpoints.Bin_Mags.get(**wildcards).output[0]
	return expand(os.path.join(checkpoint_output, "{mag}.manifest.txt"), mag=glob_wildcards(os.path.join(checkpoint_output, "{mag}.manifest.txt")).mag)

rule concatenate_mags:
	input:
		aggregate_mags
	output:
		"{pwd_directory}/mags_done.txt"
	shell:
		"""
		touch {output}
		"""

rule clean_up_files:
	input:
		mag="{pwd_directory}/mags_done.txt",
		lin="{pwd_directory}/linear_done.txt",
		circ="{pwd_directory}/circular_done.txt",
		superg="{pwd_directory}/bins_SSU_done.txt",
		bins="{pwd_directory}/bins_done.txt",
		asm="{pwd_directory}/assembly_done.txt",
		busco="{pwd_directory}/buscoAssembly.done.txt"
	output:
		"{pwd_directory}/Wolbachia.done.txt"
	shell:
		"""
		rm {input.mag} {input.lin} {input.circ} {input.superg} {input.bins} {input.busco} {input.asm}
		touch {output}
		"""
