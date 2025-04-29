rule aeron_genome_alignment:
	input:
		fastq=get_raw_fastq
	output:
		mapping_gam='{project}/{sample}/alignment/{sample}_aeron_genome_4aeron.gam',
    params:
        gfa='{project}/resources/genome.gfa',
        tmp_dir='{project}/{sample}/fusion/aeron/tmp/seedcache',
        ext_opt=config['params']['aeron_genome_alignment']
	benchmark:
        'benchmarks/{project}/aeron_genome_alignment/{sample}.txt'
    log:
        'logs/{project}/aeron_genome_alignment/{sample}.log'
	threads: config['threads']['aeron']
	shell:
		"GraphAligner -g {params.gfa} -f {input.fastq} "
        "-t {threads} "
        "--seeds-mxm-cache-prefix {params.tmp_dir} "
        "-a {output.mapping_gam}  2> {log}"



rule aeron_postprocess:
	input:
		mapping_gam='{project}/{sample}/alignment/{sample}_aeron_genome_4aeron.gam',
		fastq = get_raw_fastq
	output:
		selected_alns = '{project}/{sample}/alignment/{sample}_aeron_genome_4aeron_selected.gam',
		full_len_alns = '{project}/{sample}/alignment/{sample}_aeron_genome_4aeron_fulllength.gam',
		summary = '{project}/{sample}/fusion/aeron/{sample}_aeron_postprocess_summary.txt',
    log:
        'logs/{project}/aeron_postprocess/{sample}.log'
	benchmark:
        'benchmarks/{project}/aeron_postprocess/{sample}.txt'
    threads: config['threads']['aeron']
	shell:
		"Postprocess {input.mapping_gam} {input.fastq} {output.selected_alns} {output.full_len_alns} {output.summary} 2>{log}"


# configfile: "config.yaml"
# SCRIPTPATH=config["scripts"]
# ALIGNERBINPATH = config["binaries"]
# READFILE_FULLNAME = config["reads"]
# TRANSCRIPTFILE_FULLNAME = config["transcripts"]
# GRAPHFILE_FULLNAME = config["graph"]
# FUSION_MAX_ERROR_RATE = config["fusion_max_error_rate"]
# FUSION_MIN_SCORE_DIFFERENCE = config["fusion_min_score_difference"]
# SEEDSIZE = config["seedsize"]
# MAXSEEDHITS = config["maxseeds"]
# ALIGNERBANDWIDTH = config["aligner_bandwidth"]

# if isinstance(READFILE_FULLNAME, str): READFILE_FULLNAME = [READFILE_FULLNAME]

# READFILE_NAME = [path.split('.')[0] for path in READFILE_FULLNAME]
# TRANSCRIPTFILE_NAME = TRANSCRIPTFILE_FULLNAME.split('.')[0]
# GRAPHFILE_NAME = GRAPHFILE_FULLNAME.split('.')[0]

# def readfile_withending(wildcards):
# 	for r in READFILE_FULLNAME:
# 		if r.split('.')[0] == wildcards.reads:
# 			return "input/" + r
# 	if wildcards.reads == TRANSCRIPTFILE_NAME: return "input/" + TRANSCRIPTFILE_FULLNAME
# 	assert False

# wildcard_constraints:
# 	reads = READFILE_NAME[0],
# 	transcripts = TRANSCRIPTFILE_NAME,
# 	graph = GRAPHFILE_NAME

# rule all:
# 	input:
# 		expand("fusionoutput/reads_tofusions_onlyfusion_{reads}_{transcript}_{graph}.bam", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
# 		expand("fusionoutput/reads_tofusions_onlyfusion_{reads}_{transcript}_{graph}.bam.bai", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
# 		expand("fusionoutput/reads_tofusions_{reads}_{transcript}_{graph}.bam", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
# 		expand("fusionoutput/reads_tofusions_{reads}_{transcript}_{graph}.bam.bai", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
# 		expand("fusionoutput/fusionread_table_{reads}_{transcript}_{graph}.txt", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
# 		expand("fusionoutput/fusion_transcripts_{reads}_{transcript}_{graph}.fa", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
# 		expand("fusionoutput/ref_and_fusion_{reads}_{transcript}_{graph}.fa", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME)


# rule partial_pairs:
# 	input:
# 		readaln = "output/aln_{reads}_{graph}_secondary.gam",
# 		reads = readfile_withending
# 	output: temp("fusiontmp/pairs_{reads}_{graph}.gam")
# 	shell: "{ALIGNERBINPATH}/PickAdjacentAlnPairs {input.readaln} 20 {input.reads} {output} 100"

# rule pair_assignments:
# 	input:
# 		transcriptaln = "output/aln_{transcripts}_{graph}_full_length.gam",
# 		pairs = rules.partial_pairs.output,
# 		reads = readfile_withending
# 	output: temp("fusiontmp/matrix_{reads}_{transcripts}_{graph}.txt")
# 	shell: "{ALIGNERBINPATH}/AlignmentSubsequenceIdentity {input.transcriptaln} {input.pairs} {input.reads} 1 > {output}"

# rule exact_pair_assignments:
# 	input: "fusiontmp/matrix_{reads}_{transcripts}_{graph}.txt"
# 	output: temp("fusiontmp/exactmatrix_{reads}_{transcripts}_{graph}.txt")
# 	shell: "awk -F '\\t' '{{if ($3 == 1) print;}}' < {input} > {output}"

# rule parse_matrix:
# 	input: "fusiontmp/exactmatrix_{reads}_{transcripts}_{graph}.txt"
# 	output: temp("fusiontmp/exactparsematrix_{reads}_{transcripts}_{graph}.txt")
# 	shell: "{SCRIPTPATH}/parse_pair_matrix.py < {input} > {output}"

# rule loose_fusions:
# 	input: "fusiontmp/exactparsematrix_{reads}_{transcripts}_{graph}.txt"
# 	output: temp("fusiontmp/loose_gene_fusion_{reads}_{transcripts}_{graph}.txt")
# 	shell: "{SCRIPTPATH}/pairmatrix_get_genes.py < {input} > {output}"

# rule fusionfinder:
# 	input:
# 		graph = "input/{graph}.gfa",
# 		loose_fusions = rules.loose_fusions.output,
# 		pairmatrix = rules.exact_pair_assignments.output,
# 		transcriptaln = "output/aln_{transcripts}_{graph}_full_length.gam",
# 		reads = readfile_withending
# 	output:
# 		fusions = "fusiontmp/unfiltered_fusions_{reads}_{transcripts}_{graph}.txt",
# 		corrected = "fusiontmp/unfiltered_corrected_{reads}_{transcripts}_{graph}.txt"
# 	log:
# 		stderr = "fusiontmp/fusionfinder_stderr_{reads}_{transcripts}_{graph}.txt",
# 		stdout = "fusiontmp/fusionfinder_stdout_{reads}_{transcripts}_{graph}.txt"
# 	threads: 40
# 	shell: "/usr/bin/time -v {ALIGNERBINPATH}/FusionFinder {input.graph} {input.loose_fusions} {input.pairmatrix} {input.transcriptaln} {input.reads} 1 1.0 1 1 {threads} {output.fusions} {output.corrected} 1> {log.stdout} 2> {log.stderr}"

# rule filter_fusions:
# 	input: rules.fusionfinder.output.fusions
# 	output: "fusionoutput/fusionread_table_{reads}_{transcripts}_{graph}.txt"
# 	shell: "awk -F '\t' '{{if ($2 < {FUSION_MAX_ERROR_RATE} && $3 < -{FUSION_MIN_SCORE_DIFFERENCE}) print;}}' < {input} > {output}"

# rule fusion_transcripts:
# 	input:
# 		fusions = rules.filter_fusions.output,
# 		corrected = rules.fusionfinder.output.corrected
# 	output: "fusionoutput/fusion_transcripts_{reads}_{transcripts}_{graph}.fa"
# 	shell: "{SCRIPTPATH}/pick_fusion_exemplar.py {input.fusions} {input.corrected} > {output}"

# rule merge_ref_and_fusions:
# 	input:
# 		ref = "input/" + TRANSCRIPTFILE_FULLNAME,
# 		fusions = rules.fusion_transcripts.output
# 	output: "fusionoutput/ref_and_fusion_{reads}_{transcripts}_{graph}.fa"
# 	shell: "cat {input.ref} {input.fusions} > {output}"

# rule align_reads_to_fusions:
# 	input:
# 		refplusfusion = rules.merge_ref_and_fusions.output,
# 		reads = readfile_withending
# 	output: temp("fusiontmp/reads_tofusions_{reads}_{transcripts}_{graph}.sam")
# 	log:
# 		stderr = "fusiontmp/minimap2_stderr_{reads}_{transcripts}_{graph}.txt"
# 	threads: 40
# 	shell: "/usr/bin/time -v minimap2 --secondary=no -t {threads} -a -x map-ont {input.refplusfusion} {input.reads} 1> {output} 2> {log.stderr}"

# rule sam_to_bam:
# 	input: "fusiontmp/{filename}.sam"
# 	output: "fusionoutput/{filename}.bam"
# 	shell: "samtools view -b < {input} | samtools sort > {output}"

# rule index_bam:
# 	input: "fusionoutput/{filename}.bam"
# 	output: "fusionoutput/{filename}.bam.bai"
# 	shell: "samtools index -b {input}"

# rule fusion_support_sam:
# 	input:
# 		fusiontranscripts = rules.fusion_transcripts.output,
# 		minimapalns = rules.align_reads_to_fusions.output
# 	output:
# 		supportfile = "fusionoutput/fusion_support_{reads}_{transcripts}_{graph}.txt",
# 		fusionalnfile = temp("fusiontmp/reads_tofusions_onlyfusion_noheader_{reads}_{transcripts}_{graph}.sam")
# 	shell: "{SCRIPTPATH}/fusion_support_sam.py {input.fusiontranscripts} {input.minimapalns} 150 {output.supportfile} {output.fusionalnfile}"

# rule reheader_sam:
# 	input:
# 		fusionaln = rules.fusion_support_sam.output.fusionalnfile,
# 		header = rules.align_reads_to_fusions.output
# 	output:
# 		temp("fusiontmp/reads_tofusions_onlyfusion_{reads}_{transcripts}_{graph}.sam")
# 	run:
# 		shell("samtools view -H < {input.header} > {output}"),
# 		shell("cat {input.fusionaln} >> {output}")
















# configfile: "config.yaml"
# VGPATH = config["vgpath"]
# SCRIPTPATH=config["scripts"]
# ALIGNERBINPATH = config["binaries"]
# READFILE_FULLNAME = config["reads"]
# TRANSCRIPTFILE_FULLNAME = config["transcripts"]
# GRAPHFILE_FULLNAME = config["graph"]
# SEEDSIZE = config["seedsize"]
# MAXSEEDHITS = config["maxseeds"]
# ALIGNERBANDWIDTH = config["aligner_bandwidth"]
# GTFPATH = config["gtffile"]
# ALIGNMENTSELECTION = config["alignment_selection"]
# ECUTOFF = config["alignment_E_cutoff"]

# #input files at top: check them!

# # all input files must be in the folder ./input/
# # use the full file name, including file ending

# # input splice graph
# # format must be .vg
# graph: HumanUpdated38V5.gfa

# # reference transcripts
# # format can be either fasta/fastq, gzipped or not
# transcripts: HumanTranscriptsReduced.fa
# # sequenced reads
# # format can be either fasta/fastq, gzipped or not
# # for more files, add them in new lines starting with "- "
# # NOTE: the file names without ending must be unique! You cannot have eg. reads.fq and reads.fa
# reads: humanNA12878dRNA.fasta

# #optional params below: default values will probably work

# #size of the seed hits. Fewer means more accurate but slower alignments.
# seedsize: 17
# #max number of seeds. Fewer means faster but more inaccurate alignment
# maxseeds: 15

# fusion_max_error_rate: 0.2
# fusion_min_score_difference: 200

# #options: --greedy-length --greedy-score --greedy-E --schedule-score --schedule-length --schedule-inverse-E-product --schedule-inverse-E-sum
# alignment_selection: 
# alignment_E_cutoff: 1

# #bandwidth for the aligner. Higher means more accurate but slower alignment.
# aligner_bandwidth: 35

# gtffile: HomoSapiens.gtf

# #file paths

# # https://bitbucket.org/dilipdurai/aeron/
# scripts: /MMCI/MS/LongRNA/work/mikko/AERON_runs/run_20200515_threeprime/Aeron/AeronScripts
# # https://github.com/maickrau/GraphAligner
# binaries: /MMCI/MS/LongRNA/work/mikko/AERON_runs/run_20200515_threeprime/Aeron/Binaries
# # needed to convert mummer seeds to .gam seeds
# vgpath: /MMCI/MS/LongRNA/work/mikko/vg/vg


# if isinstance(READFILE_FULLNAME, str): READFILE_FULLNAME = [READFILE_FULLNAME]

# READFILE_NAME = [path.split('.')[0] for path in READFILE_FULLNAME]
# TRANSCRIPTFILE_NAME = TRANSCRIPTFILE_FULLNAME.split('.')[0]
# GRAPHFILE_NAME = GRAPHFILE_FULLNAME.split('.')[0]




# configfile: "config.yaml"
# SCRIPTPATH=config["scripts"]
# ALIGNERBINPATH = config["binaries"]
# READFILE_FULLNAME = config["reads"]
# TRANSCRIPTFILE_FULLNAME = config["transcripts"]
# GRAPHFILE_FULLNAME = config["graph"]
# FUSION_MAX_ERROR_RATE = config["fusion_max_error_rate"]
# FUSION_MIN_SCORE_DIFFERENCE = config["fusion_min_score_difference"]
# SEEDSIZE = config["seedsize"]
# MAXSEEDHITS = config["maxseeds"]
# ALIGNERBANDWIDTH = config["aligner_bandwidth"]

# if isinstance(READFILE_FULLNAME, str): READFILE_FULLNAME = [READFILE_FULLNAME]

# READFILE_NAME = [path.split('.')[0] for path in READFILE_FULLNAME]
# TRANSCRIPTFILE_NAME = TRANSCRIPTFILE_FULLNAME.split('.')[0]
# GRAPHFILE_NAME = GRAPHFILE_FULLNAME.split('.')[0]

# def readfile_withending(wildcards):
# 	for r in READFILE_FULLNAME:
# 		if r.split('.')[0] == wildcards.reads:
# 			return "input/" + r
# 	if wildcards.reads == TRANSCRIPTFILE_NAME: return "input/" + TRANSCRIPTFILE_FULLNAME
# 	assert False

# wildcard_constraints:
# 	reads = READFILE_NAME[0],
# 	transcripts = TRANSCRIPTFILE_NAME,
# 	graph = GRAPHFILE_NAME

# rule all:
# 	input:
# 		expand("fusionoutput/reads_tofusions_onlyfusion_{reads}_{transcript}_{graph}.bam", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
# 		expand("fusionoutput/reads_tofusions_onlyfusion_{reads}_{transcript}_{graph}.bam.bai", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
# 		expand("fusionoutput/reads_tofusions_{reads}_{transcript}_{graph}.bam", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
# 		expand("fusionoutput/reads_tofusions_{reads}_{transcript}_{graph}.bam.bai", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
# 		expand("fusionoutput/fusionread_table_{reads}_{transcript}_{graph}.txt", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
# 		expand("fusionoutput/fusion_transcripts_{reads}_{transcript}_{graph}.fa", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
# 		expand("fusionoutput/ref_and_fusion_{reads}_{transcript}_{graph}.fa", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME)

# rule align_with_secondaries:
# 	input:
# 		graph = "input/{graph}.gfa",
# 		reads = readfile_withending
# 	output:
# 		"output/aln_{reads}_{graph}_secondary.gam"
# 	log:
# 		stdout = "tmp/aligner_stdout_{reads}_{graph}.txt",
# 		stderr = "tmp/aligner_stderr_{reads}_{graph}.txt"
# 	threads: 40
# 	shell:
# 		"/usr/bin/time -v {ALIGNERBINPATH}/Aligner --all-alignments -g {input.graph} -f {input.reads} --try-all-seeds --seeds-mxm-length {SEEDSIZE} --seeds-mem-count {MAXSEEDHITS} --seeds-mxm-cache-prefix tmp/seeds_{wildcards.graph}_index -a {output} -t {threads} -b {ALIGNERBANDWIDTH} --E-cutoff 1 1> {log.stdout} 2> {log.stderr}"

# rule partial_pairs:
# 	input:
# 		readaln = "output/aln_{reads}_{graph}_secondary.gam",
# 		reads = readfile_withending
# 	output: temp("fusiontmp/pairs_{reads}_{graph}.gam")
# 	shell: "{ALIGNERBINPATH}/PickAdjacentAlnPairs {input.readaln} 20 {input.reads} {output} 100"

# rule pair_assignments:
# 	input:
# 		transcriptaln = "output/aln_{transcripts}_{graph}_full_length.gam",
# 		pairs = rules.partial_pairs.output,
# 		reads = readfile_withending
# 	output: temp("fusiontmp/matrix_{reads}_{transcripts}_{graph}.txt")
# 	shell: "{ALIGNERBINPATH}/AlignmentSubsequenceIdentity {input.transcriptaln} {input.pairs} {input.reads} 1 > {output}"

# rule exact_pair_assignments:
# 	input: "fusiontmp/matrix_{reads}_{transcripts}_{graph}.txt"
# 	output: temp("fusiontmp/exactmatrix_{reads}_{transcripts}_{graph}.txt")
# 	shell: "awk -F '\\t' '{{if ($3 == 1) print;}}' < {input} > {output}"

# rule parse_matrix:
# 	input: "fusiontmp/exactmatrix_{reads}_{transcripts}_{graph}.txt"
# 	output: temp("fusiontmp/exactparsematrix_{reads}_{transcripts}_{graph}.txt")
# 	shell: "{SCRIPTPATH}/parse_pair_matrix.py < {input} > {output}"

# rule loose_fusions:
# 	input: "fusiontmp/exactparsematrix_{reads}_{transcripts}_{graph}.txt"
# 	output: temp("fusiontmp/loose_gene_fusion_{reads}_{transcripts}_{graph}.txt")
# 	shell: "{SCRIPTPATH}/pairmatrix_get_genes.py < {input} > {output}"

# rule fusionfinder:
# 	input:
# 		graph = "input/{graph}.gfa",
# 		loose_fusions = rules.loose_fusions.output,
# 		pairmatrix = rules.exact_pair_assignments.output,
# 		transcriptaln = "output/aln_{transcripts}_{graph}_full_length.gam",
# 		reads = readfile_withending
# 	output:
# 		fusions = "fusiontmp/unfiltered_fusions_{reads}_{transcripts}_{graph}.txt",
# 		corrected = "fusiontmp/unfiltered_corrected_{reads}_{transcripts}_{graph}.txt"
# 	log:
# 		stderr = "fusiontmp/fusionfinder_stderr_{reads}_{transcripts}_{graph}.txt",
# 		stdout = "fusiontmp/fusionfinder_stdout_{reads}_{transcripts}_{graph}.txt"
# 	threads: 40
# 	shell: "/usr/bin/time -v {ALIGNERBINPATH}/FusionFinder {input.graph} {input.loose_fusions} {input.pairmatrix} {input.transcriptaln} {input.reads} 1 1.0 1 1 {threads} {output.fusions} {output.corrected} 1> {log.stdout} 2> {log.stderr}"

# rule filter_fusions:
# 	input: rules.fusionfinder.output.fusions
# 	output: "fusionoutput/fusionread_table_{reads}_{transcripts}_{graph}.txt"
# 	shell: "awk -F '\t' '{{if ($2 < {FUSION_MAX_ERROR_RATE} && $3 < -{FUSION_MIN_SCORE_DIFFERENCE}) print;}}' < {input} > {output}"

# rule fusion_transcripts:
# 	input:
# 		fusions = rules.filter_fusions.output,
# 		corrected = rules.fusionfinder.output.corrected
# 	output: "fusionoutput/fusion_transcripts_{reads}_{transcripts}_{graph}.fa"
# 	shell: "{SCRIPTPATH}/pick_fusion_exemplar.py {input.fusions} {input.corrected} > {output}"

# rule merge_ref_and_fusions:
# 	input:
# 		ref = "input/" + TRANSCRIPTFILE_FULLNAME,
# 		fusions = rules.fusion_transcripts.output
# 	output: "fusionoutput/ref_and_fusion_{reads}_{transcripts}_{graph}.fa"
# 	shell: "cat {input.ref} {input.fusions} > {output}"

# rule align_reads_to_fusions:
# 	input:
# 		refplusfusion = rules.merge_ref_and_fusions.output,
# 		reads = readfile_withending
# 	output: temp("fusiontmp/reads_tofusions_{reads}_{transcripts}_{graph}.sam")
# 	log:
# 		stderr = "fusiontmp/minimap2_stderr_{reads}_{transcripts}_{graph}.txt"
# 	threads: 40
# 	shell: "/usr/bin/time -v minimap2 --secondary=no -t {threads} -a -x map-ont {input.refplusfusion} {input.reads} 1> {output} 2> {log.stderr}"

# rule sam_to_bam:
# 	input: "fusiontmp/{filename}.sam"
# 	output: "fusionoutput/{filename}.bam"
# 	shell: "samtools view -b < {input} | samtools sort > {output}"

# rule index_bam:
# 	input: "fusionoutput/{filename}.bam"
# 	output: "fusionoutput/{filename}.bam.bai"
# 	shell: "samtools index -b {input}"

# rule fusion_support_sam:
# 	input:
# 		fusiontranscripts = rules.fusion_transcripts.output,
# 		minimapalns = rules.align_reads_to_fusions.output
# 	output:
# 		supportfile = "fusionoutput/fusion_support_{reads}_{transcripts}_{graph}.txt",
# 		fusionalnfile = temp("fusiontmp/reads_tofusions_onlyfusion_noheader_{reads}_{transcripts}_{graph}.sam")
# 	shell: "{SCRIPTPATH}/fusion_support_sam.py {input.fusiontranscripts} {input.minimapalns} 150 {output.supportfile} {output.fusionalnfile}"

# rule reheader_sam:
# 	input:
# 		fusionaln = rules.fusion_support_sam.output.fusionalnfile,
# 		header = rules.align_reads_to_fusions.output
# 	output:
# 		temp("fusiontmp/reads_tofusions_onlyfusion_{reads}_{transcripts}_{graph}.sam")
# 	run:
# 		shell("samtools view -H < {input.header} > {output}"),
# 		shell("cat {input.fusionaln} >> {output}")

