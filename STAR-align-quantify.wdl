workflow SRA_Downloads{

    File RSEM_index
    File genomeGTF
    String title
    String links
    Int read_length
    String library_type
    File in_fastq
    
    if(library_type=="SINGLE"){
       	call align_fastq_single{
         input:
            fastq = in_fastq,
            genomeGTF = genomeGTF,
            title = title,
            index_files = RSEM_index,
            links=links
	} 
    	call bam_output_single{
        input:
            AlignedSAM = align_fastq_single.AlignedSAM,
            title =  title

    }
    	output{
    
    	File? ReadsPerGene_s = align_fastq_single.ReadsPerGene
        File? SortedBAM_s = bam_output_single.BAM
        File? IndexBAM_s = bam_output_single.BAMind
    		
    }
    }
    if(library_type=="PAIRED"){
       	call align_fastq_paired{
         input:
            fastq = in_fastq,
            genomeGTF = genomeGTF,
            title = title,
            index_files = RSEM_index,
            links=links
	}   call bam_output_paired{
         input:
            AlignedSAM = align_fastq_paired.AlignedSAM,
            title =  title

    }
    	output{
    
    	File? ReadsPerGene_p = align_fastq_paired.ReadsPerGene
        File? SortedBAM_p = bam_output_paired.BAM
        File? IndexBAM_p = bam_output_paired.BAMind
    		
    }
    }
    
}

 task align_fastq_paired {

    String title
    File fastq
    File index_files
    File genomeGTF
    String links
        	
    
    command <<<
    
    echo ${links} | tr " " "\n" > gs_links.txt
    for value in `cat gs_links.txt`; do echo $value | cut -d/ -f5; done > gs_title.txt
    paste gs_links.txt  gs_title.txt > gs_out.txt
    SRRID=`cat gs_out.txt | awk '{print $2}' | uniq`
    echo $SRRID
    
    # Unpack index and fastq files 
    
    tar -xvf ${index_files}
    tar -xvf ${fastq}
    # Align fastq
    
    STAR --runThreadN 4 \
    	--runMode alignReads \
        --outFilterMultimapNmax 1 \
        --outFilterMatchNmin 35 \
        --quantMode GeneCounts \
        --twopassMode Basic \
        --outFileNamePrefix ${title} \
        --sjdbGTFfile ${genomeGTF} \
        --genomeDir out/GenomeDir \
        --readFilesIn paired_merged_fastq/$SRRID_1.fastq paired_merged_fastq/$SRRID_2.fastq
    
    
    head ${title}ReadsPerGene.out.tab
    
    
        >>>
    output {
    
    	File ReadsPerGene = "${title}ReadsPerGene.out.tab"
        File AlignedSAM = "${title}Aligned.out.sam"
    
    }

    runtime {
        disks: 			"local-disk 125 HDD"
        cpu: 			16
        memory: 		"72 GB"
        docker: 		"gargk/sra"
        preemptible: 	3
    }
    
 }
 task align_fastq_single {

    String title
    File fastq
    File index_files
    File genomeGTF
    String links
            	
    
    command <<<
    
    # get SRR ID
    echo ${links} | tr " " "\n" > gs_links.txt
    for value in `cat gs_links.txt`; do echo $value | cut -d/ -f5; done > gs_title.txt
    paste gs_links.txt  gs_title.txt > gs_out.txt
    SRRID=`cat gs_out.txt | awk '{print $2}' | uniq`
    echo $SRRID
    # Unpack index and fastq files 
    
    tar -xvf ${index_files}
    tar -xvf ${fastq}
    # Align fastq
    
    STAR --runThreadN 4 \
    	--runMode alignReads \
        --outFilterMultimapNmax 1 \
        --outFilterMatchNmin 35 \
        --quantMode GeneCounts \
        --twopassMode Basic \
        --outFileNamePrefix ${title} \
        --sjdbGTFfile ${genomeGTF} \
        --genomeDir out/GenomeDir \
        --readFilesIn single_merged_fastq/$${title}.fastq
    
    
    head ${title}ReadsPerGene.out.tab
    
    
        >>>
    output {
    
    	File ReadsPerGene = "${title}ReadsPerGene.out.tab"
        File AlignedSAM = "${title}Aligned.out.sam"
    
    }

    runtime {
        disks: 			"local-disk 125 HDD"
        cpu: 			16
        memory: 		"72 GB"
        docker: 		"gargk/sra"
        preemptible: 	3
    }
    
 }


 task bam_output_single {

    File AlignedSAM
    String title
            	
    
    command <<<
    
    	samtools view -b -S ${AlignedSAM} > ${title}.bam
		samtools sort ${title}.bam -o st.${title}.bam
        samtools index st.${title}.bam 
        

        >>>
    output {
    
    	File BAM = "st.${title}.bam"
        File BAMind = "st.${title}.bam.bai"
    	
    }

    runtime {
        disks: 			"local-disk 100 HDD"
        cpu: 			4
        memory: 		"16 GB"
        docker: 		"jweinstk/samtools"
        preemptible: 	3
    }
    
 }
 task bam_output_paired {

    File AlignedSAM
    String title
            	
    
    command <<<
    
    	samtools view -b -S ${AlignedSAM} > ${title}.bam
		samtools sort ${title}.bam -o st.${title}.bam
        samtools index st.${title}.bam 
        

        >>>
    output {
    
    	File BAM = "st.${title}.bam"
        File BAMind = "st.${title}.bam.bai"
    	
    }

    runtime {
        disks: 			"local-disk 100 HDD"
        cpu: 			4
        memory: 		"16 GB"
        docker: 		"jweinstk/samtools"
        preemptible: 	3
    }
    
 }