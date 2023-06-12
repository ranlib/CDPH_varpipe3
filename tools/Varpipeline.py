#! /usr/bin/env python
"""
determine variants
"""
import subprocess
import os
from datetime import datetime
import yaml


class snp:
    """get variants"""

    def __init__(self, input, outdir, reference, reference_name, name, paired, input2, verbose, threads, argString):
        i = datetime.now()
        self.name = name  # sample name
        self.flog = "Output_" + i.strftime("%m_%d_%Y")
        self.qlog = os.path.join(self.flog, "QC")
        self.fOut = os.path.join(self.flog, outdir)
        self.input = input
        self.outdir = os.path.join(self.fOut, "tmp")
        self.tmp = os.path.join(self.outdir, "tmp")
        self.clockwork = os.path.join(self.fOut, "clockwork")
        self.trimmomatic = os.path.join(self.fOut, "trimmomatic")
        self.paired = paired
        self.input2 = input2
        self.verbose = verbose
        self.reference = reference
        self.reference_name = reference_name
        self.__lineage = os.path.join(self.fOut, self.name + ".lineage_report.txt")
        self.__finalVCF = ""
        self.__fullVCF = ""
        self.__annotation = ""
        self.__full_annotation = ""
        self.__inter_annotation = ""
        self.__final_annotation = ""
        self.__mixed = ""
        self._low = ""
        self.__exception = ""
        self.__threads = threads

        # Create the output directory, and start the log file.
        self.__logged = False
        if not os.path.isfile(self.flog):
            self.__CallCommand("mkdir", ["mkdir", self.flog])
        if not os.path.isfile(self.qlog):
            self.__CallCommand("mkdir", ["mkdir", self.qlog])

        self.__CallCommand("mkdir", ["mkdir", self.fOut])

        self.__CallCommand("mkdir", ["mkdir", "-p", self.tmp])
        self.__CallCommand("mkdir", ["mkdir", "-p", self.trimmomatic])
        self.__CallCommand("mkdir", ["mkdir", "-p", self.clockwork])
        self.__log = os.path.join(self.fOut, self.name + ".log")
        self.__qlog = os.path.join(self.qlog, "QC.log")

        cwd = os.getcwd()
        configname = os.path.join(os.path.dirname(__file__), "config.yml")
        fh10 = open(configname, "r")
        fh20 = open(os.path.join(cwd, "config.yml"), "w")
        for lines in fh10:
            line = lines.rstrip("\r\n")
            if line.startswith("tools") or line.startswith("scripts") or line.startswith("other"):
                print(line, file=fh20)
            elif "threads" in line:
                print(line, file=fh20)
            else:
                lined = line.split(":")
                print(lined[0] + ":" + " " + os.path.dirname(__file__) + lined[1][1:], file=fh20)
        fh10.close()
        fh20.close()

        with open(os.path.join(cwd, "config.yml"), "r") as ymlfile:
            cfg = yaml.load(ymlfile, Loader=yaml.Loader)
        self.__logFH = open(self.__log, "w")
        self.__logFH2 = open(self.__qlog, "a")
        self.__logFH.write(argString + "\n\n")
        self.__logged = True

        # Format Validation
        self.__nextflow = cfg["tools"]["nextflow"]
        self.__remove_contam = cfg["tools"]["remove_contam"]
        self.__ref_fasta = cfg["tools"]["ref_fasta"]
        self.__ref_metadata = cfg["tools"]["ref_metadata"]
        self.__trimmomatic = cfg["tools"]["trimmomatic"]
        # Mapping
        self.__bwa = cfg["tools"]["bwa"]
        self.__samtools = cfg["tools"]["samtools"]
        # Picard-Tools
        self.__picard = cfg["tools"]["picard"]
        # SNP / InDel Calling
        self.__gatk = cfg["tools"]["gatk"]
        # Other
        self.__bcftools = cfg["tools"]["bcftools"]
        self.__bedtools = cfg["tools"]["bedtools"]
        self.__annotator = cfg["tools"]["annotator"]
        self.__parser = cfg["scripts"]["parser"]
        self.__creater = cfg["scripts"]["creater"]
        self.__interpreter = cfg["scripts"]["interpreter"]
        self.__included = cfg["scripts"]["included"]
        self.__reported = cfg["scripts"]["reported"]
        self.__stats_estimator = cfg["scripts"]["stats_estimator"]
        self.__genome_stats_estimator = cfg["scripts"]["genome_stats_estimator"]
        self.__bedlist_amp = cfg["scripts"]["bedlist_amp"]
        self.__bedlist_one = cfg["scripts"]["bedlist_one"]
        self.__bedlist_two = cfg["scripts"]["bedlist_two"]
        self.__bedstruct = cfg["scripts"]["bedstruct"]
        self.__snpeff_database = cfg["scripts"]["snpeff_database"]
        self.__target_estimator = cfg["scripts"]["target_coverage_estimator"]
        self.__genome_coverage_estimator = cfg["scripts"]["genome_coverage_estimator"]
        self.mutationloci = cfg["scripts"]["mutationloci"]
        self.__lineage_parser = cfg["scripts"]["lineage_parser"]
        self.__lineages = cfg["scripts"]["lineages"]
        self.__structparser = cfg["scripts"]["struct_parser"]
        self.__create_report = cfg["scripts"]["create_report"]
        self.__print_report = cfg["scripts"]["print_report"]
        self.__threads = cfg["other"]["threads"]

    def __CallCommand(self, program, command):
        """Allows execution of a simple command."""
        out = ""
        err = ""
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()

        if isinstance(program, list):
            o = open(program[1], "wt")
            o.write(out.decode())
            o.close()
            out = ""
            program = program[0]

        if self.__logged:
            self.__logFH.write("---[ " + program + " ]---\n")
            self.__logFH.write("Command: \n" + " ".join(command) + "\n\n")
            if out:
                self.__logFH.write("Standard Output: \n" + out.decode() + "\n\n")
            if err:
                self.__logFH.write("Standard Error: \n" + err.decode() + "\n\n")
        return 1

    def runClockwork(self):
        """Clockwork Decontamination"""
        self.__ifVerbose("Performing clockwork decontamination.")
        if self.paired:
            self.__CallCommand("nextflow remove contamination", [self.__nextflow, "run", self.__remove_contam, "--ref_fasta", self.__ref_fasta, "--ref_metadata_tsv", self.__ref_metadata, "--reads_in1", self.input, "--reads_in2", self.input2, "--outprefix", self.clockwork + "/" + self.name, "--mapping_threads", self.__threads])
            self.input = self.clockwork + "/" + self.name + ".remove_contam.1.fq.gz"
            self.input2 = self.clockwork + "/" + self.name + ".remove_contam.2.fq.gz"

    def runTrimmomatic(self):
        """QC Trimmomatic"""
        self.__ifVerbose("Performing trimmomatic trimming.")
        trimlog = self.trimmomatic + "/trimLog.txt"
        if self.paired:
            paired1 = os.path.join(self.trimmomatic, self.name + "_paired_1.fastq.gz")
            paired2 = os.path.join(self.trimmomatic, self.name + "_paired_2.fastq.gz")
            unpaired1 = os.path.join(self.trimmomatic, self.name + "_unpaired_1.fastq.gz")
            unpaired2 = os.path.join(self.trimmomatic, self.name + "_unpaired_2.fastq.gz")
            command = f"trimmomatic PE -threads {self.__threads} -trimlog {trimlog} {self.input} {self.input2} {paired1} {unpaired1} {paired2} {unpaired2} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40"
            self.__CallCommand("trimmomatic", command.split())
            # self.__CallCommand("rm", ["rm", self.trimmomatic + "/" + self.name + "_unpaired_1.fastq.gz", self.trimmomatic + "/" + self.name + "_unpaired_2.fastq.gz"])
            self.input = paired1
            self.input2 = paired2
        else:
            trimmed = os.path.join(self.trimmomatic, self.name + "_paired.fastq.gz")
            command = f"trimmomatic SE -threads {self.__threads} -trimlog {trimlog} {self.input} {trimmed} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40"
            self.__CallCommand("trimmomatic", command.split())
            self.input = trimmed

    def runBWA(self, bwa):
        """Align reads against the reference using bwa."""
        self.__ranBWA = True
        self.__ifVerbose("Running BWA.")
        self.__logFH.write("########## Running BWA. ##########\n")
        bwaOut = os.path.join(self.outdir, "bwa")
        self.__CallCommand("mkdir", ["mkdir", "-p", bwaOut])
        self.__ifVerbose("   Building BWA index.")
        self.__bwaIndex(bwaOut + "/index")
        self.__alnSam = bwaOut + "/bwa.sam"
        self.__bwaLongReads(bwaOut)
        self.__ifVerbose("   Process Alignments")
        self.__processAlignment()

    def __bwaIndex(self, out):
        """Make an index of the given reference genome."""
        self.__CallCommand("mkdir", ["mkdir", "-p", out])
        self.__CallCommand("cp", ["cp", self.reference, out + "/ref.fa"])
        self.reference = out + "/ref.fa"
        self.__CallCommand("bwa index", ["bwa", "index", self.reference])
        if os.path.exists(out + "/ref.dict"):
            os.remove(out + "/ref.dict")
        self.__CallCommand("CreateSequenceDictionary", ["gatk", "CreateSequenceDictionary", "-R", self.reference, "-O", out + "/ref.dict"])
        self.__CallCommand("samtools faidx", ["samtools", "faidx", self.reference])

    def __bwaLongReads(self, out):
        """Make use of bwa mem"""
        read_group = "@RG\\tID:" + self.name + "\\tSM:" + self.name + "\\tPL:ILLUMINA"
        if self.paired:
            self.__ifVerbose("   Running BWA mem on paired end reads.")
            command = f"bwa mem -t {str(self.__threads)} -K 10000000 -c 100 -M -T 50 -R {read_group} -o {self.__alnSam} {self.reference} {self.input} {self.input2}"
            self.__CallCommand("bwa mem", command.split())
        else:
            self.__ifVerbose("   Running BWA mem on single end reads.")
            command = f"bwa mem -t {str(self.__threads)} -K 10000000 -c 100 -M -T 50 -R {read_group} -o {self.__alnSam} {self.reference} {self.input}"
            self.__CallCommand("bwa mem", command.split())

    def __processAlignment(self):
        """Filter alignment using GATK and Picard-Tools"""
        self.__ifVerbose("Filtering alignment with GATK and Picard-Tools.")
        self.__logFH.write("########## Filtering alignment with GATK and Picard-Tools. ##########\n")
        GATKdir = os.path.join(self.outdir, "GATK")
        self.__CallCommand("mkdir", ["mkdir", "-p", GATKdir])
        samDir = os.path.join(self.outdir, "SamTools")
        self.__CallCommand("mkdir", ["mkdir", "-p", samDir])

        """ Convert SAM to BAM"""
        if self.__ranBWA:
            self.__ifVerbose("   Running SamFormatConverter.")
            self.__CallCommand(
                "SamFormatConverter",
                [
                    "gatk",
                    "SamFormatConverter",
                    "-INPUT",
                    self.__alnSam,
                    "-VALIDATION_STRINGENCY",
                    "LENIENT",
                    "-OUTPUT",
                    GATKdir + "/GATK.bam",
                ],
            )
        else:
            self.__CallCommand("cp", ["cp", self.__alnSam, GATKdir + "/GATK.bam"])

        """ Run mapping Report and Mark duplicates using Picard-Tools"""
        self.__ifVerbose("   Running SortSam.")
        self.__CallCommand("SortSam", ["gatk", "SortSam", "-INPUT", GATKdir + "/GATK.bam", "-SORT_ORDER", "coordinate", "-OUTPUT", GATKdir + "/GATK_s.bam", "-VALIDATION_STRINGENCY", "LENIENT", "-TMP_DIR", self.tmp])
        self.__ifVerbose("   Running MarkDuplicates.")
        self.__CallCommand("MarkDuplicates", ["gatk", "MarkDuplicates", "-INPUT", GATKdir + "/GATK_s.bam", "-OUTPUT", GATKdir + "/GATK_sdr.bam", "-METRICS_FILE", GATKdir + "/MarkDupes.metrics", "-ASSUME_SORTED", "true", "-REMOVE_DUPLICATES", "false", "-VALIDATION_STRINGENCY", "LENIENT"])
        self.__ifVerbose("   Running BuildBamIndex.")
        self.__CallCommand("BuildBamIndex", ["gatk", "BuildBamIndex", "-INPUT", GATKdir + "/GATK_sdr.bam", "-VALIDATION_STRINGENCY", "LENIENT"])
        self.__CallCommand("samtools view", ["samtools", "view", "-c", "-o", samDir + "/unmapped.txt", GATKdir + "/GATK_sdr.bam"])

        """ Filter out unmapped reads """
        self.__finalBam = os.path.join(self.fOut, self.name + "_sdrcsm.bam")
        self.__ifVerbose("   Running samtools view.")
        self.__CallCommand("samtools view", ["samtools", "view", "-bhF", "4", "-o", self.__finalBam, GATKdir + "/GATK_sdr.bam"])
        self.__ifVerbose("   Running BuildBamIndex.")
        self.__CallCommand("BuildBamIndex", ["gatk", "BuildBamIndex", "-INPUT", self.__finalBam, "-VALIDATION_STRINGENCY", "LENIENT"])
        self.__ifVerbose("   Running samtools view")
        self.__CallCommand("rm", ["rm", "-r", self.tmp])
        self.__CallCommand("samtools view", ["samtools", "view", "-c", "-o", samDir + "/mapped.txt", self.__finalBam])

    def runCoverage(self):
        """Run genome Coverage Statistics"""

        self.__ifVerbose("Running target Coverage Statistics")
        samDir = self.outdir + "/SamTools"
        i = datetime.now()

        self.__CallCommand("samtools depth", ["samtools", "depth", "-o", samDir + "/coverage.txt", "-a", self.__finalBam])
        self.__CallCommand(["bedtools coverage", samDir + "/bed_amp_coverage.txt"], ["bedtools", "coverage", "-b", self.__finalBam, "-a", self.__bedlist_amp])
        self.__CallCommand(["sort", samDir + "/bed_amp_sorted_coverage.txt"], ["sort", "-nk", "6", samDir + "/bed_amp_coverage.txt"])
        self.__CallCommand(["bedtools coverage", samDir + "/bed_1_coverage.txt"], ["bedtools", "coverage", "-b", self.__finalBam, "-a", self.__bedlist_one])
        self.__CallCommand(["bedtools coverage", samDir + "/bed_2_coverage.txt"], ["bedtools", "coverage", "-b", self.__finalBam, "-a", self.__bedlist_two])
        self.__CallCommand(["sort", samDir + "/bed_1_sorted_coverage.txt"], ["sort", "-nk", "6", samDir + "/bed_1_coverage.txt"])
        self.__CallCommand(["sort", samDir + "/bed_2_sorted_coverage.txt"], ["sort", "-nk", "6", samDir + "/bed_2_coverage.txt"])
        self.__CallCommand(["target region coverage estimator", samDir + "/target_region_coverage_amp.txt"], [self.__target_estimator, self.__bedlist_amp, samDir + "/coverage.txt", self.name])
        self.__CallCommand(["sort", self.fOut + "/" + self.name + "_target_region_coverage.txt"], ["sort", "-nk", "3", samDir + "/target_region_coverage_amp.txt"])
        self.__CallCommand(["genome stats estimator", samDir + "/" + self.name + "_genome_stats.txt"], [self.__genome_stats_estimator, samDir + "/coverage.txt", self.name])
        self.__CallCommand(["genome coverage estimator", samDir + "/genome_region_coverage_1.txt"], [self.__genome_coverage_estimator, samDir + "/bed_1_sorted_coverage.txt", samDir + "/coverage.txt", self.name])
        self.__CallCommand(["genome coverage estimator", samDir + "/genome_region_coverage_2.txt"], [self.__genome_coverage_estimator, samDir + "/bed_2_sorted_coverage.txt", samDir + "/coverage.txt", self.name])
        self.__CallCommand(["cat", samDir + "/genome_region_coverage.txt"], ["cat", samDir + "/genome_region_coverage_1.txt", samDir + "/genome_region_coverage_2.txt"])
        self.__CallCommand(["sort", self.fOut + "/" + self.name + "_genome_region_coverage.txt"], ["sort", "-nk", "3", samDir + "/genome_region_coverage.txt"])
        self.__CallCommand("sed", ["sed", "-i", "1d", self.fOut + "/" + self.name + "_genome_region_coverage.txt"])
        self.__CallCommand(["structural variant detector", self.fOut + "/" + self.name + "_structural_variants.txt"], [self.__structparser, self.__bedstruct, self.fOut + "/" + self.name + "_genome_region_coverage.txt", samDir + "/coverage.txt", self.name])
        self.__CallCommand(["stats estimator", self.fOut + "/" + self.name + "_stats.txt"], [self.__stats_estimator, samDir + "/unmapped.txt", samDir + "/mapped.txt", self.fOut + "/" + self.name + "_target_region_coverage.txt", self.name, samDir + "/" + self.name + "_genome_stats.txt"])
        statsOut = self.fOut + "/" + self.name + "_stats.txt"
        fh20 = open(statsOut, "r")
        for lines in fh20:
            if lines.startswith("Sample ID"):
                continue
            fields = lines.rstrip("\r\n").split("\t")
            if float(fields[2]) < 90.0 or float(fields[3]) < 30 or float(fields[4]) < 90.0:
                self.__logFH2.write(i.strftime("%Y/%m/%d %H:%M:%S") + "\t" + "Sample:" + "\t" + self.name + "\t" + "failed QC checks\n")
                # self.__CallCommand('rm', ['rm', '-r', self.outdir])
                # self.__CallCommand('rm', ['rm',  self.__finalBam])
                # self.__CallCommand('rm', ['rm', '-r', self.trimmomatic])
                # self.__CallCommand('rm', ['rm', '-r', self.clockwork])
                # self.__CallCommand('mv', ['mv', self.fOut, self.qlog])
                # sys.exit(1)

    def runGATK(self):
        """Variant Callers"""
        if os.path.isfile(self.__finalBam):
            self.__ifVerbose("Calling SNPs/InDels with GATK.")
            self.__logFH.write("########## Calling SNPs/InDels with Mutect2. ##########\n")
            GATKdir = os.path.join(self.outdir, "GATK")
            samDir = os.path.join(self.outdir, "SamTools")

            """ Set final VCF file. """
            if not self.__finalVCF:
                self.__finalVCF = os.path.join(GATKdir, self.name + "_filter.vcf")

            if not self.__fullVCF:
                self.__fullVCF = os.path.join(GATKdir, self.name + "_full_filter.vcf")

            """ Call SNPs/InDels with Mutect2 """
            self.__ifVerbose("   Running Mutect2.")
            self.__CallCommand("Mutect2", ["gatk", "Mutect2", "-R", self.reference, "-I", self.__finalBam, "-O", GATKdir + "/mutect.vcf", "--max-mnp-distance", "2", "-L", self.__included])

            self.__CallCommand("Mutect2", ["gatk", "Mutect2", "-R", self.reference, "-I", self.__finalBam, "-O", GATKdir + "/full_mutect.vcf", "--max-mnp-distance", "2"])

            self.__CallCommand("LeftAlignAndTrimVariants", ["gatk", "LeftAlignAndTrimVariants", "-R", self.reference, "-V", GATKdir + "/mutect.vcf", "-O", GATKdir + "/gatk_mutect.vcf", "--split-multi-allelics"])
            self.__CallCommand("mv", ["mv", GATKdir + "/gatk_mutect.vcf", GATKdir + "/mutect.vcf"])

            self.__CallCommand("LeftAlignAndTrimVariants", ["gatk", "LeftAlignAndTrimVariants", "-R", self.reference, "-V", GATKdir + "/full_mutect.vcf", "-O", GATKdir + "/full_gatk_mutect.vcf", "--split-multi-allelics"])
            self.__CallCommand("mv", ["mv", GATKdir + "/full_gatk_mutect.vcf", GATKdir + "/full_mutect.vcf"])

            self.__CallCommand("FilterMutectCalls", ["gatk", "FilterMutectCalls", "-R", self.reference, "-V", GATKdir + "/mutect.vcf", "--min-reads-per-strand", "1", "--min-median-read-position", "10", "--min-allele-fraction", "0.01", "--microbial-mode", "true", "-O", self.__finalVCF])

            self.__CallCommand("FilterMutectCalls", ["gatk", "FilterMutectCalls", "-R", self.reference, "-V", GATKdir + "/full_mutect.vcf", "--min-reads-per-strand", "1", "--min-median-read-position", "10", "--min-allele-fraction", "0.01", "--microbial-mode", "true", "-O", self.__fullVCF])

        else:
            print("<E> no bam file found")

    def annotateVCF(self):
        """Annotate the final VCF file"""
        if self.__finalVCF:
            self.__ifVerbose("Annotating final VCF.")
            self.__annotation = os.path.join(self.fOut, self.name + "_DR_loci_raw_annotation.vcf")
            self.__CallCommand(["SnpEff", self.__annotation], ["snpEff", "-nodownload", "-noLog", "-noStats", "-c", self.__snpeff_database, self.reference_name, self.__finalVCF])

            self.__ifVerbose("Parsing final Annotation.")
            # translate vcf -> tsv file
            self.__CallCommand(["create annotation", self.fOut + "/" + self.name + "_DR_loci_annotation.txt"], [self.__creater, self.__annotation, self.name])
            # vcf -> tsv + some manipulation
            self.__CallCommand(["parse annotation", self.fOut + "/" + self.name + "_DR_loci_Final_annotation.txt"], [self.__parser, self.__annotation, self.mutationloci, self.name])

        if self.__fullVCF:
            self.__ifVerbose("Annotating full VCF.")
            self.__full_annotation = os.path.join(self.fOut, self.name + "_full_raw_annotation.vcf")
            self.__CallCommand(["SnpEff", self.__full_annotation], ["snpEff", "-nodownload", "-noLog", "-noStats", "-c", self.__snpeff_database, self.reference_name, self.__fullVCF])

            self.__ifVerbose("Parsing full Annotation.")
            self.__CallCommand(["create annotation", self.fOut + "/" + self.name + "_full_annotation.txt"], [self.__creater, self.__full_annotation, self.name])
            self.__CallCommand(["parse annotation", self.fOut + "/" + self.name + "_full_Final_annotation.txt"], [self.__parser, self.__full_annotation, self.mutationloci, self.name])
        else:
            self.__ifVerbose("Use SamTools, GATK, or Freebayes to annotate the final VCF.")
        # cwd = os.getcwd()
        # self.__CallCommand('rm', ['rm',  cwd + "/snpEff_genes.txt"])
        # self.__CallCommand('rm', ['rm',  cwd + "/snpEff_summary.html"])

    def runLineage(self):
        """Run lineage Analysis"""
        self.__ifVerbose("Running Lineage Analysis")
        self.__full_final_annotation = os.path.join(self.fOut, self.name + "_full_Final_annotation.txt")
        self.__CallCommand(["lineage parsing", self.fOut + "/" + self.name + "_Lineage.txt"], [self.__lineage_parser, self.__lineages, self.__full_final_annotation, self.__lineage, self.name])

    def runPrint(self):
        """Print analysis report"""
        self.__ifVerbose("Printing report")
        self.__CallCommand(["create summary report", self.fOut + "/" + self.name + "_summary.txt"], [self.__create_report, self.fOut + "/" + self.name + "_stats.txt", self.fOut + "/" + self.name + "_target_region_coverage.txt", self.fOut + "/" + self.name + "_DR_loci_Final_annotation.txt"])
        self.__CallCommand(["run interpretation report", self.fOut + "/" + self.name + "_interpretation.txt"], [self.__interpreter, self.__reported, self.fOut + "/" + self.name + "_summary.txt", self.fOut + "/" + self.name + "_structural_variants.txt", self.fOut + "/" + self.name + "_DR_loci_annotation.txt", self.fOut + "/" + self.name + "_target_region_coverage.txt", self.name])
        self.__CallCommand("print pdf report", [self.__print_report, self.fOut + "/" + self.name + "_summary.txt", self.fOut + "/" + self.name + "_report.pdf"])

    def cleanUp(self):
        """Clean up the temporary files, and move them to a proper folder."""
        i = datetime.now()
        cwd = os.getcwd()
        self.__CallCommand("rm", ["rm", "-r", self.outdir])
        self.__CallCommand("rm", ["rm", self.fOut + "/" + self.name + "_sdrcsm.bai"])
        self.__CallCommand("rm", ["rm", self.__finalBam])
        self.__CallCommand("rm", ["rm", "-r", self.trimmomatic])
        self.__CallCommand("rm", ["rm", "-r", self.clockwork])
        self.__CallCommand("rm", ["rm", os.path.join(cwd, "config.yml")])
        self.__CallCommand("rm", ["rm", "-r", os.path.join(cwd, "work")])
        if os.path.isfile(self.fOut + "/" + self.name + ".log"):
            self.__logFH.close()
            self.__logged = False
        fh4 = open(self.__log, "r")
        for line in fh4:
            lines = line.rstrip("\r\n")
            if "Exception" in lines:
                self.__exception = "positive"
        fh4.close()
        if self.__exception == "positive":
            self.__logFH2.write(i.strftime("%Y/%m/%d %H:%M:%S") + "\t" + "Input:" + "\t" + self.name + "\t" + "Exception in analysis\n")
            self.__CallCommand("mv", ["mv", self.fOut, self.qlog])

        if os.path.isfile(self.fOut + "/" + self.name + "_stats.txt"):
            if os.path.getsize(self.fOut + "/" + self.name + "_stats.txt") < 1:
                self.__logFH2.write(i.strftime("%Y/%m/%d %H:%M:%S") + "\t" + "Input:" + "\t" + self.name + "\t" + "Exception in analysis\n")
                self.__CallCommand("mv", ["mv", self.fOut, self.qlog])

    def __ifVerbose(self, msg):
        """If verbose print a given message"""
        if self.verbose:
            print(msg)
