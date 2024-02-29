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

rule all:
	input:
		expand('{pwd_dir}/symbiont_presence.txt', pwd_dir=config["pwd_directory"]),
		expand('{pwd_dir}/assembly_done.txt', pwd_dir=config["pwd_directory"]),
		expand('{pwd_dir}/buscoAssembly.done.txt', pwd_dir=config["pwd_directory"]),
		expand('{pwd_dir}/busco.done.txt', pwd_dir=config["pwd_directory"]),
		expand('{pwd_dir}/bins_done.txt', pwd_dir=config["pwd_directory"]),
		expand('{pwd_dir}/total.done.txt', pwd_dir=config["pwd_directory"])

rule CheckPresence:
	"""
	check if this sample contains Wolbachia or Spiroplasma
	"""
	input:
		taxfile = expand("{tax}", tax=config["taxfile"])
	output:
		wolpresence = expand("{pwd_directory}/symbiont_presence.txt", pwd_directory=config["pwd_directory"])
	shell:
			"""
			if grep -q 'Wolbachia;' {input.taxfile}; then
				echo "Wolbachia,Anaplasmataceae" >> {output.wolpresence}
			fi

			if grep -q 'Mesenet;' {input.taxfile}; then
				echo "Mesenet,Anaplasmataceae" >> {output.wolpresence}
			fi

			if grep -q 'Spiroplasma;' {input.taxfile}; then
				echo "Spiroplasma,Spiroplasmataceae" >> {output.wolpresence}
			fi

			if grep -q 'Cardinium;' {input.taxfile}; then
				echo "Cardinium,Amoebophilaceae" >> {output.wolpresence}
			fi

			if grep -q 'Rickettsiella;' {input.taxfile}; then
				echo "Cardinium,Coxiellaceae" >> {output.wolpresence}
			fi	

			if grep -q 'Arsenophonus;' {input.taxfile}; then
				echo "Arsenophonus,Morganellaceae" >> {output.wolpresence}
			fi

			if grep -q 'Lariskella;' {input.taxfile}; then
				echo "Lariskella,Midichloriaceae" >> {output.wolpresence}
			fi

			if grep -q 'Sulcia;' {input.taxfile}; then
				echo "Sulcia,Flavobacteriales" >> {output.wolpresence}
			fi

			if grep -q 'Tremblaya;' {input.taxfile}; then
				echo "Tremblaya,Betaproteobacteria" >> {output.wolpresence}
			fi	

			if grep -q 'Zinderia;' {input.taxfile}; then
				echo "Zinderia,Oxalobacteraceae" >> {output.wolpresence}
			fi	

			if grep -q 'Buchnera;' {input.taxfile}; then
				echo "Buchnera,Erwiniaceae" >> {output.wolpresence}
			fi	

			if grep -q 'Portiera;' {input.taxfile}; then
				echo "Portiera,Halomonadaceae" >> {output.wolpresence}
			fi	

			if grep -q 'Carsonella;' {input.taxfile}; then
				echo "Carsonella,Halomonadaceae" >> {output.wolpresence}
			fi																	
			"""

checkpoint CheckPresenceCobionts:
	"""
	check if this sample contains Wolbachia or Spiroplasma
	"""
	input:
		wolpresence = expand("{pwd_directory}/symbiont_presence.txt", pwd_directory=config["pwd_directory"])
	output:
		famdir = directory(expand("{pwd_directory}/families/", pwd_directory=config["pwd_directory"])),
	shell:
			"""
			if [ ! -d {output.famdir} ]; then
				mkdir {output.famdir}
			fi

			if [ -s {input.wolpresence} ]; then
				while read p
            	do
            		shortname=`echo $p | cut -d, -f1`
					family=`echo $p | cut -d, -f2`       
                	echo $p > {output.famdir}/genus.$family.txt
            	done < {input.wolpresence}
			fi
			"""

rule RunHifiasm:
	"""
	check if hifiasm was already run, if not re-try
	"""  
	input:
		wolpresence = expand("{pwd_directory}/symbiont_presence.txt", pwd_directory=config["pwd_directory"]),
		fasta_reads = expand("{pwd}", pwd=config["workingdirectory"]),
		#generafiles = glob.glob("{pwd_directory}/families/genus.*.txt")
		generafiles = "{pwd_directory}/families/genus.{family}.txt"
	params:
		assemblyprefix = "{pwd_directory}/{family}/hifiasm/hifiasm",
		dirname = directory(expand("{pwd}",pwd=config["workingdirectory"])),
		fam = "{family}"
	output:
		completed = "{pwd_directory}/{family}.assembly_done.txt",
		dirname2 = directory("{pwd_directory}/{family}/hifiasm"),
		gfa = "{pwd_directory}/{family}/hifiasm/hifiasm.p_ctg.gfa",
		fasta = "{pwd_directory}/{family}/hifiasm/hifiasm.p_ctg.fasta",
		fai = "{pwd_directory}/{family}/hifiasm/hifiasm.p_ctg.fasta.fai"
	threads: 10
	conda: "envs/hifiasm_seqtk.yaml"
	shell:
            """
			if [ ! -d {output.dirname2} ]; then
				mkdir {output.dirname2}
				cp {input.fasta_reads}/{params.fam}/kraken.fa {pwd_dir}/{params.fam}/hifiasm/reads.fa
			fi
            if [ -s {input.fasta_reads}/{params.fam}/kraken.fa ] && [ ! -s {params.dirname}/{params.fam}/hifiasm/hifiasm.p_ctg.noseq.gfa ] ; then
				echo {input.fasta_reads}/{params.fam}/kraken.fa
				echo {params.dirname}/{params.fam}/hifiasm/hifiasm.p_ctg.noseq.gfa
				linecount=$(grep -c '>' < {input.fasta_reads}/{params.fam}/kraken.fa)
				if [ $linecount -ge 50000 ]; then
					seqtk sample {input.fasta_reads}/{params.fam}/kraken.fa 50000 > {pwd_dir}/{params.fam}/hifiasm/reads.fa
				fi
				hifiasm -o {params.assemblyprefix} -t {threads} {pwd_dir}/{params.fam}/hifiasm/reads.fa -D 10 -l 1 -s 0.999 || true
				mv {pwd_dir}/{params.fam}/hifiasm/hifiasm.bp.p_ctg.gfa  {output.gfa}
				awk '/^S/{{print ">"$2"\\n"$3}}' {output.gfa} | fold > {output.fasta} || true
				faidx {output.fasta}
			else
				cp {params.dirname}/{params.fam}/hifiasm/hifiasm.p_ctg.gfa {output.gfa} 
				cp {params.dirname}/{params.fam}/hifiasm/hifiasm.p_ctg.fasta {output.fasta} 
				cp {params.dirname}/{params.fam}/hifiasm/hifiasm.p_ctg.fasta.fai {output.fai}
			fi
			touch {output.completed}
			"""

def aggregate_asms(wildcards):
	checkpoint_output=checkpoints.CheckPresenceCobionts.get(**wildcards).output[0]
	return expand(os.path.join("{pwd_dir}/{family}.assembly_done.txt"), pwd_dir=config["pwd_directory"],family=glob_wildcards(os.path.join(checkpoint_output, "genus.{family}.txt")).family )

rule concatenate_asms:
	input:
		aggregate_asms
	output:
		"{pwd_directory}/assembly_done.txt"
	shell:
		"""
		touch {output}
		"""

rule RunBuscoAssembly:
    """
	Detect number of BUSCO genes per contig
	"""
	input:
		hifiasm_done = "{pwd_directory}/{family}.assembly_done.txt",
		circgenome = "{pwd_directory}/{family}/hifiasm/hifiasm.p_ctg.fasta",
		generafiles = "{pwd_directory}/families/genus.{family}.txt"
	params:
		buscodir = directory(expand("{pwd}/",pwd=config["workingdirectory"])),
		buscodir2 = "{pwd_directory}/{family}/buscoAssembly",
		buscodbs = "{pwd_directory}/{family}/info_dbs_assembly.txt",
		buscoini = "{pwd_directory}/{family}/config_busco_assembly.ini",
		fam = "{family}"
	output:
		completed = "{pwd_directory}/{family}.buscoAssembly.done.txt",
		summary = "{pwd_directory}/{family}/buscoAssembly/completeness_per_contig.txt"
	conda: "envs/busco.yaml"
	threads: 10
	shell:
            """
			if [ ! -d {params.buscodir2} ]; then
				mkdir {params.buscodir2}
			fi

			if [ -f {params.buscodir}/{params.fam}/buscoAssembly/busco/run*/full_table.tsv ] || [ -f {params.buscodir}/{params.fam}/buscoAssembly/full_table.tsv ]; then
				if [ -f {params.buscodir}/{params.fam}/buscoAssembly/full_table.tsv ]; then
					cp {params.buscodir}/{params.fam}/buscoAssembly/full_table.tsv {params.buscodir2}
				else
					cp {params.buscodir}/{params.fam}/buscoAssembly/busco/run*/full_table.tsv {params.buscodir2}
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

def aggregate_buscos(wildcards):
	checkpoint_output=checkpoints.CheckPresenceCobionts.get(**wildcards).output[0]
	return expand(os.path.join("{pwd_dir}/{family}.buscoAssembly.done.txt"), pwd_dir=config["pwd_directory"],family=glob_wildcards(os.path.join(checkpoint_output, "genus.{family}.txt")).family )

rule concatenate_buscos:
	input:
		aggregate_buscos
	output:
		"{pwd_directory}/buscoAssembly.done.txt"
	shell:
		"""
		touch {output}
		"""

rule RunBusco:
    """
	Detect number of BUSCO genes per contig
	"""
	input:
		hifiasm_done = "{pwd_directory}/{family}.assembly_done.txt",
		busco_done = "{pwd_directory}/{family}.buscoAssembly.done.txt",
		generafiles = "{pwd_directory}/families/genus.{family}.txt"
	params:
		circgenome =  "{pwd_directory}/{family}/hifiasm/hifiasm.p_ctg.fasta",
		buscodir = directory(expand("{pwd}",pwd=config["workingdirectory"])),
		buscodir2 = directory("{pwd_directory}/{family}/busco"),
		buscodbs = "{pwd_directory}/{family}/info_dbs_assembly.txt",
		buscoini = "{pwd_directory}/{family}/config_busco.ini",
		fam = "{family}"
	output:
		completed = "{pwd_directory}/{family}.busco.done.txt",
		summary = "{pwd_directory}/{family}/busco/completeness_per_contig.txt"
	conda: "envs/busco.yaml"
	threads: 10
	shell:
            """
			if [ ! -d {params.buscodir2} ]; then
				mkdir {params.buscodir2}
			fi

			if [ -f {params.buscodir}/{params.fam}/busco/busco/run*/full_table.tsv ] || [ -f {params.buscodir}/{params.fam}/busco/full_table.tsv ]; then
				if [ -f {params.buscodir}/{params.fam}/busco/full_table.tsv ]; then
					cp {params.buscodir}/{params.fam}/busco/full_table.tsv {params.buscodir2}
				else
					cp {params.buscodir}/{params.fam}/busco/busco/run*/full_table.tsv {params.buscodir2}
				fi
			else
				if [ -s {pwd}/{params.fam}.finalassembly.fa ] ; then
					busco --list-datasets > {params.buscodbs}
                	python {scriptdir}/BuscoConfig.py -f {pwd}/{params.fam}.finalassembly.fa -d {params.buscodir2} -dl busco_data/ -c {threads} -db {params.buscodbs} -o {params.buscoini}
                	busco --config {params.buscoini} -f || true
					mv {params.buscodir2}/busco/run*/full_table.tsv {params.buscodir2} 
				else
					touch {params.buscodir2}/full_table.tsv
				fi
			fi
            touch {output.completed}
			if [ -s {pwd}/{params.fam}/{params.fam}.finalassembly.fa ] ; then
				python {scriptdir}/ParseBuscoTableMapping.py -d {params.buscodir2} -i {pwd}/{params.fam}/{params.fam}.finalassembly.fa -o {output.summary}
			elif [ -s {pwd}/{params.fam}/{params.fam}.ctgs.fa ] ; then
				python {scriptdir}/ParseBuscoTableMapping.py -d {params.buscodir2} -i {pwd}/{params.fam}/{params.fam}.ctgs.fa -o {output.summary}
			else
				touch {output.summary} 
			fi
            """

def aggregate_buscos2(wildcards):
	checkpoint_output=checkpoints.CheckPresenceCobionts.get(**wildcards).output[0]
	return expand(os.path.join("{pwd_dir}/{family}.busco.done.txt"), pwd_dir=config["pwd_directory"],family=glob_wildcards(os.path.join(checkpoint_output, "genus.{family}.txt")).family )

rule concatenate_buscos2:
	input:
		aggregate_buscos2
	output:
		"{pwd_directory}/busco.done.txt"
	shell:
		"""
		touch {output}
		"""

rule WolbachiaQuality:
	"""
	check quality of assemblies, create metadata file and chrom list file (if necessary). if quality not sufficient, create mag files
	"""
	input:
		generafiles = "{pwd_directory}/families/genus.{family}.txt",
		hifiasm_done = "{pwd_directory}/assembly_done.txt",
		busco_asm_done = "{pwd_directory}/buscoAssembly.done.txt",
		busco_done = "{pwd_directory}/busco.done.txt",
		summary = "{pwd_directory}/{family}/buscoAssembly/completeness_per_contig.txt",
		orig_dir = expand("{pwd}",pwd=config["workingdirectory"])
	output:
		completed = "{pwd_directory}/{family}.bins_done.txt",
		binlist = "{pwd_directory}/{family}/bin_list.txt"
	shell:
		"""
		python {scriptdir}/WolbachiaQuality.py -d {input.summary} -t {tolid} -l {output.binlist} -g {gfa} -d2 {input.orig_dir} -o {outdir} -f {input.generafiles}
		touch {output.completed}
		"""

checkpoint Bins:
	"""
	Move all fasta bins to directory supergroup
	"""
	input:
		bins_quality = "{pwd_directory}/{family}.bins_done.txt",
		binlist = "{pwd_directory}/{family}/bin_list.txt"
	params:
		fam = "{family}"
	output:
		superdir = directory("{pwd_directory}/{family}/supergroup/")
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
				mv {pwd_dir}/{params.fam}/$shortname.fa {output.superdir}/bin.$shortname.fa
			fi
        done < {input.binlist}
		"""

def aggregate_bins(wildcards):
	#checkpoint_output=checkpoints.Bins.get(**wildcards).output[0]
	checkpoint_output=checkpoints.CheckPresenceCobionts.get(**wildcards).output[0]
	#return expand(os.path.join(checkpoint_output, "{family}.bins_done.txt"), family=glob_wildcards(os.path.join(pwd_dir, "genus.{family}.txt")).family)
	return expand(os.path.join("{pwd_dir}/{family}.bins_done.txt"), pwd_dir=config["pwd_directory"],family=glob_wildcards(os.path.join(checkpoint_output, "genus.{family}.txt")).family )

rule concatenate_bins:
	input:
		aggregate_bins
	output:
		"{pwd_directory}/bins_done.txt"
	shell:
		"""
		touch {output}
		"""

rule HMMscan_SSU:
	"""
	Run HMMscan with prokaryotic+viral HMM (RF00177+RF01959)
	"""
	input:
		fasta = "{pwd_directory}/{family}/supergroup/bin.{bin}.fa"
	output:
		domfile = "{pwd_directory}/{family}/supergroup/{bin}.ProkSSU.domout",
		log = "{pwd_directory}/{family}/supergroup/{bin}.HMMscan.log"
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
		dom = "{pwd_directory}/{family}/supergroup/{bin}.ProkSSU.domout", 
		fasta = "{pwd_directory}/{family}/supergroup/bin.{bin}.fa"
	output:
		readsinfo = "{pwd_directory}/{family}/supergroup/{bin}.ProkSSU.readsinfo",
		fasta16SLoci = "{pwd_directory}/{family}/supergroup/{bin}.ProkSSU.fa",
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
		fasta16SLoci = "{pwd_directory}/{family}/supergroup/{bin}.ProkSSU.fa"
	output:
		aln16SLoci = "{pwd_directory}/{family}/supergroup/{bin}.ProkSSU.aln.fa",
		treefile = "{pwd_directory}/{family}/supergroup/{bin}.ProkSSU.treefile"
	params:
		fam = "{family}"
	conda: "envs/treebuild.yaml"
	threads: 10
	shell:
		"""
		if [ -s {input.fasta16SLoci} ] && [ {params.fam} = 'Anaplasmataceae' ]; then
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
                donefile = temporary("{pwd_dir}/taxdownload.done.txt")
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
		treefile = "{pwd_directory}/{family}/supergroup/{bin}.ProkSSU.treefile",
		fasta16SLoci = "{pwd_directory}/{family}/supergroup/{bin}.ProkSSU.fa",
		taxfile = expand("{tax}", tax=config["taxfile"]),
		taxdone = "{pwd_directory}/taxdownload.done.txt",
		generafiles = "{pwd_directory}/families/genus.{family}.txt"
	params:
		binname = "{bin}",
		fam = "{family}"
	output:
		supergroup_name =  "{pwd_directory}/{family}/supergroup/{bin}.SSU.supergroup.txt"
	conda: "envs/ete3.yaml"
	shell:	
		"""
		if [ -s {input.fasta16SLoci} ] ; then
			python {scriptdir}/WolbachiaSupergroup.py -t {input.treefile} -c {input.fasta16SLoci} -f {input.generafiles} -b {params.binname} -s {host} -o {output.supergroup_name} -sts {stsfile} -i {tolid} -na {ncbi_taxdir}/names.dmp -no {ncbi_taxdir}/nodes.dmp -ta {input.taxfile}
			speciesname=` cat {output.supergroup_name} | cut -f1`
			today="$(date +'%Y%m%d')"
			if [ -s {outdir}/{params.binname}."$today"/{params.binname}.yaml ]; then
				echo "species: "$speciesname >>  {outdir}/{params.binname}."$today"/{params.binname}.yaml
			fi
			#cp {output.supergroup_name} {outdir}/{params.binname}.speciesname.txt
			#write code here to add novel name to file
			cat {output.supergroup_name} >> "{pwd_dir}/{params.fam}/taxid.nr.txt"
		else
			touch {output.supergroup_name} "{pwd_dir}/{params.fam}/taxid.nr.txt"
		fi
		"""

def aggregate_supergroups(wildcards):
	checkpoint_output=checkpoints.Bins.get(**wildcards).output[0]
	#return expand("{workingdirectory}/{family}/supergroup/{bin}.SSU.supergroup.txt", workingdirectory=config["pwd_directory"], bin=glob_wildcards(os.path.join(checkpoint_output, 'bin.{bin}.fa')).bin, family=glob_wildcards(os.path.join(checkpoint_output2, "genus.{family}.txt")).family )
	return expand(os.path.join(checkpoint_output, "{bin}.SSU.supergroup.txt"), bin=glob_wildcards(os.path.join(checkpoint_output, 'bin.{bin}.fa')).bin)
	#return expand("{workingdirectory}/{family}/supergroup/{bin}.SSU.supergroup.txt", workingdirectory=config["pwd_directory"], bin=glob_wildcards(os.path.join(checkpoint_output, 'bin.{bin}.fa')).bin, family=glob_wildcards(os.path.join({workingdirectory}, "families/genus.{family}.txt")).family )


checkpoint concatenate_supergroups:
	input:
		aggregate_supergroups
	output:
		reqfile="{pwd_directory}/{family}.taxid_request.novel.txt",
		donefile="{pwd_directory}/{family}/bins_SSU_done.txt"
	params:
		fam = "{family}"
	shell:
		"""
		#write code here to submit email to ENA taxid request
		cut -f1-5 {pwd_dir}/{params.fam}/taxid.nr.txt | sort | uniq > {pwd_dir}/{params.fam}/taxid_request.txt

		#add script double-check it doesn't exist
		python {scriptdir}/taxid_requests.py -i {pwd_dir}/{params.fam}/taxid_request.txt -o {output.reqfile}

		rm {pwd_dir}/{params.fam}/taxid.nr.txt
		sh {scriptdir}/email_request_taxids.sh {output.reqfile} {email}
		touch {output.donefile}
		"""

checkpoint Bin_Circular:
	"""
	Move all fasta bins to directory supergroup
	"""
	input:
		binsdone = "{pwd_directory}/{family}/bins_SSU_done.txt",
		bins_quality = "{pwd_directory}/{family}/bins_done.txt",
		binlist = "{pwd_directory}/{family}/bin_list.txt"
	output:
		superdir = directory("{pwd_directory}/{family}/circular/")
	params:
		fam = "{family}"
	shell:
		"""
		if [ ! -d {output.superdir} ]; then
			mkdir {output.superdir}
        fi
 		while read p
		do
			shortname=`echo $p | cut -d, -f1`
			echo {pwd_dir}/{params.fam}/$shortname.list
			if [ -f {pwd_dir}/{params.fam}/$shortname.list ] && grep "Circular" {pwd_dir}/{params.fam}/$shortname.list ; then
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
		fasta = "{pwd_directory}/{family}/supergroup/bin.{circ}.fa",
		binname= "{pwd_directory}/{family}/circular/bin.{circ}.txt"
	output:
		circ_fa = "{pwd_directory}/{family}/circular/{circ}.fa"
	params:
		circname = "{circ}",
		fam = "{family}"
	threads: 10
	conda: "envs/prokka.yaml"
	shell:
		"""
		prokka --cpus {threads} --outdir {pwd_dir}/{params.fam}/circular/{params.circname}.prokka --prefix {params.circname} {input.fasta}
		python {scriptdir}/RotateRevComp.py -gff {pwd_dir}/{params.fam}/circular/{params.circname}.prokka/{params.circname}.gff -fa {input.fasta} -o {output.circ_fa} -f {params.fam}
		today="$(date +'%Y%m%d')"
		rm {outdir}"/"{params.circname}"."$today"/"{params.circname}".fa.gz"
		cp {output.circ_fa} {outdir}"/"{params.circname}"."$today
		gzip {outdir}"/"{params.circname}"."$today"/"{params.circname}".fa"

		sh {scriptdir}/cobiont_curation_request.sh {outdir}"/"{params.circname}"."$today"/"{params.circname}".yaml" {email}
		"""

def aggregate_circs(wildcards):
	checkpoint_output=checkpoints.Bin_Circular.get(**wildcards).output[0]
	return expand(os.path.join(checkpoint_output, "{bin}.fa"), bin=glob_wildcards(os.path.join(checkpoint_output, 'bin.{bin}.txt')).bin)
	#return expand(os.path.join(checkpoint_output, "{circ}.fa"), circ=glob_wildcards(os.path.join(checkpoint_output, 'bin.{circ}.txt')).circ)

rule concatenate_circular:
	input:
		aggregate_circs
	output:
		"{pwd_directory}/{family}/circular_done.txt"
	shell:
		"""
		touch {output}
		"""

checkpoint Bin_Linear:
	"""
	Move all fasta bins to directory supergroup
	"""
	input:
		binsdone = "{pwd_directory}/{family}/bins_SSU_done.txt",
		bins_quality = "{pwd_directory}/{family}/bins_done.txt",
		binlist = "{pwd_directory}/{family}/bin_list.txt"
	output:
		superdir = directory("{pwd_directory}/{family}/linear/")
	params:
		fam = "{family}"
	shell:
		"""
		if [ ! -d {output.superdir} ]; then
			mkdir {output.superdir}
        fi
 		while read p
		do
			shortname=`echo $p | cut -d, -f1`
			today="$(date +'%Y%m%d')"
			if [ -f {pwd_dir}/{params.fam}/$shortname.list ] && grep "Linear" {pwd_dir}/{params.fam}/$shortname.list ; then
				touch {output.superdir}/bin.$shortname.txt
				sh {scriptdir}/cobiont_curation_request.sh {outdir}"/"$shortname"."$today"/"$shortname".yaml" {email}
			elif [ ! -f {pwd_dir}/{params.fam}/$shortname.list ] && [ ! -f {pwd_dir}/{params.fam}/$shortname.manifest.txt ] && [ -f {pwd_dir}/{params.fam}/$shortname.metadata.txt ] ; then
				touch {output.superdir}/bin.$shortname.txt
				sh {scriptdir}/cobiont_curation_request.sh {outdir}"/"$shortname"."$today"/"$shortname".yaml" {email}
			fi
        done < {input.binlist}
		"""

def aggregate_lins(wildcards):
	checkpoint_output=checkpoints.Bin_Linear.get(**wildcards).output[0]
	#return expand("{workingdirectory}/{family}/linear/bin.{bin}.fa", workingdirectory=config["pwd_directory"], bin=glob_wildcards(os.path.join(checkpoint_output, 'bin.{bin}.fa')).bin, family=glob_wildcards(os.path.join({workingdirectory}, "families/genus.{family}.txt")).family )
	return expand(os.path.join(checkpoint_output, "bin.{circ}.txt"), circ=glob_wildcards(os.path.join(checkpoint_output, 'bin.{circ}.txt')).circ)

rule concatenate_linear:
	input:
		aggregate_lins
	output:
		"{pwd_directory}/{family}/linear_done.txt"
	shell:
		"""
		touch {output}
		"""

checkpoint Bin_Mags:
	"""
	Move all fasta bins to directory supergroup
	"""
	input:
		binsdone = "{pwd_directory}/{family}/bins_SSU_done.txt",
		bins_quality = "{pwd_directory}/{family}/bins_done.txt",
		binlist = "{pwd_directory}/{family}/bin_list.txt",
	output:
		superdir = directory("{pwd_directory}/{family}/mag/")
	params:
		fam = "{family}"
	shell:
		"""
		if [ ! -d {output.superdir} ]; then
			mkdir {output.superdir}
        fi
 		while read p
		do
			shortname=`echo $p | cut -d, -f1`
			if [ -f {pwd_dir}/{params.fam}/$shortname.manifest.txt ] ; then
				mv {pwd_dir}/{params.fam}/$shortname.manifest.txt {output.superdir}
				#cp {output.superdir}/$shortname.manifest.txt {outdir}
				#cp {pwd_dir}/{params.fam}/supergroup/bin.$shortname.fa {outdir}/$shortname.mag.fa 
			fi
        done < {input.binlist}
		"""

def aggregate_mags(wildcards):
	checkpoint_output=checkpoints.Bin_Mags.get(**wildcards).output[0]
	return expand(os.path.join(checkpoint_output, "{mag}.manifest.txt"), mag=glob_wildcards(os.path.join(checkpoint_output, "{mag}.manifest.txt")).mag)
	#return expand("{workingdirectory}/{family}/mag/{bin}.manifest.txt", workingdirectory=config["pwd_directory"], bin=glob_wildcards(os.path.join(checkpoint_output, 'bin.{bin}.fa')).bin, family=glob_wildcards(os.path.join({workingdirectory}, "families/genus.{family}.txt")).family )

rule concatenate_mags:
	input:
		aggregate_mags
	output:
		"{pwd_directory}/{family}/mags_done.txt"
	shell:
		"""
		touch {output}
		"""

rule clean_up_files:
	input:
		mag="{pwd_directory}/{family}/mags_done.txt",
		lin="{pwd_directory}/{family}/linear_done.txt",
		circ="{pwd_directory}/{family}/circular_done.txt",
		superg="{pwd_directory}/{family}/bins_SSU_done.txt"
	output:
		"{pwd_directory}/{family}.total.done.txt"
	shell:
		"""
		rm {input.superg}
		touch {output}
		"""

def aggregate_clean_by_fam(wildcards):
	checkpoint_output=checkpoints.CheckPresenceCobionts.get(**wildcards).output[0]
	return expand(os.path.join("{pwd_dir}/{family}.total.done.txt"), pwd_dir=config["pwd_directory"], family=glob_wildcards(os.path.join(checkpoint_output, "genus.{family}.txt")).family )

rule concatenate_clean_total:
	input:
		aggregate_clean_by_fam
	output:
		"{pwd_directory}/total.done.txt"
	shell:
		"""
		rm {pwd_dir}/bins_done.txt {pwd_dir}/assembly_done.txt {pwd_dir}/buscoAssembly.done.txt
		rm {pwd_dir}/*done.txt {pwd_dir}/*/*done.txt
		touch {output}
		"""