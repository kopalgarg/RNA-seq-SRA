workflow SRA_Downloads{

    File kallisto_index
    File genomeGTF
    String title
    String links
    Int read_length
    String library_type
     
    call download_sra_fastq{
     	 input:
     	 	title = title,
            links = links
   }
   	call quantify{
         input:
            fastq = download_sra_fastq.fastq,
            genomeGTF = genomeGTF,
            title = title,
            index_files = kallisto_index,
            read_length = read_length,
            library_type = library_type
}

	output{
   
    		File counts = quantify.abundance
            File BAMIndex = quantify.pseudoalignmentsBAMIdx
            File BAM = quantify.pseudoalignmentsBAM
    }
}

 task download_sra_fastq {

    String title
    String links
            	    
    command <<<
    
    # Download SRA 
    
       wget ${links} 
       ls -lah
       cat SRR* > ${title}".sra"
       
    # Convert to FASTQ
    
      echo "fastq-dump"
      fastq-dump -O . ${title}".sra"
       
       
        >>>
    output {
    
    	File fastq = "${title}.fastq"
    
    }

    runtime {
        disks: 			"local-disk 75 HDD"
        cpu: 			4
        memory: 		"16 GB"
        docker: 		"flemoine/sratoolkit"
        preemptible: 	3
    }
    
 }

 task quantify {

    String title
    File fastq
    File index_files
    File genomeGTF
    Int read_length
    String library_type
            	
    
    command <<<
    
    # unpack index files 
    
    tar -xvf ${index_files}
    
    # quantify
        
    cd out
    
    library_type=${library_type}
    
	if [[ "$library_type" == *[sS][iI][nN][gG][lL][eE]* ]]
    then 
    	kallisto quant -i index -o ${title} -t 4 --single -l ${read_length} -s 2 --gtf ${genomeGTF} --genomebam ${fastq}
    else 
   		kallisto quant -i index -o ${title} -t 4 ${fastq}
    fi

    mv ${title}/pseudoalignments.bam.bai ../${title}_pseudoalignments.bam.bai
    mv ${title}/pseudoalignments.bam ../${title}_pseudoalignments.bam
    mv ${title}/abundance.tsv ../${title}_abundance.tsv
    
    >>>
    output {
    	
        File pseudoalignmentsBAM = "${title}_pseudoalignments.bam"
        File pseudoalignmentsBAMIdx = "${title}_pseudoalignments.bam.bai"
        File abundance = "${title}_abundance.tsv"
   
    }

    runtime {
        disks: 			"local-disk 125 HDD"
        cpu: 			16
        memory: 		"32 GB"
        docker: 		"zlskidmore/kallisto"
        preemptible: 	3
    }
    
 }


