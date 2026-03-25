###################################
####### phasing
###################################
configfile: "config/config.yaml"

## sample params
person = config["patient_id"]
CCS_Samps = config["sample_id"]
Reference_samp = config["reference_sample_id"]

## path params
software_Path = config["software_dir"]
iSNVtable_P = config["isnv_table_dir"]
refFasta_path = config["reference_fasta_dir"]
bam_path = config["bam_dir"]
OUT_path = config["output_dir"]

rule all:
    input:
        expand("{OUT_path}/readHap.{samp}_2_{Ref_samp}/2.read-haplotype.Stat.txt",
               samp=[CCS_Samps], Ref_samp=Reference_samp, OUT_path=OUT_path),

rule read_LongReadbam_file:
    input:
        reliablePosiVcf=expand("{iSNVtable_P}/{person}.reliable.locus.CCS.table.txt",
                               iSNVtable_P=iSNVtable_P, person=person),
        BamF=expand("{bam_path}/{samp}/{samp}_2_{Ref_samp}.sort.bam",
                    bam_path=bam_path, samp=[CCS_Samps], Ref_samp=Reference_samp),
        refSampGnm=expand("{refFasta_path}/{Ref_samp}.fasta",
                          refFasta_path=refFasta_path, Ref_samp=Reference_samp),
    output:
        SNV_type=expand("{OUT_path}/readHap.{samp}_2_{Ref_samp}/2.read-haplotype.txt",
                        OUT_path=OUT_path, samp=[CCS_Samps], Ref_samp=Reference_samp),
        SNV_Type_Stat=expand("{OUT_path}/readHap.{samp}_2_{Ref_samp}/2.read-haplotype.Stat.txt",
                             OUT_path=OUT_path, samp=[CCS_Samps], Ref_samp=Reference_samp),
    params:
        scriptPath=config["software_dir"],
        outP=expand("{OUT_path}/readHap.{samp}_2_{Ref_samp}",
                    OUT_path=OUT_path, samp=[CCS_Samps], Ref_samp=Reference_samp),
        ONT_PosiIndex_step="200",
    shell:
        "python {params.scriptPath}/function_read_ONTbam.py -b {input.BamF} -v {input.reliablePosiVcf} "
        "-r {input.refSampGnm} -o {params.outP} -s {params.ONT_PosiIndex_step}"
