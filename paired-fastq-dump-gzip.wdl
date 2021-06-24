workflow SRA_Downloads{
    String title
    String buckets
    call make_text{
         input:
            buckets = buckets,
            title = title    
    }
    scatter(srr in make_text.temp){
    	call srr_id {
        	input:
            	bucket = srr
        
        }
    	call sra_fastq {
    		input:
    			bucket = srr,
                srr = srr_id.id
    	}
    }
    call merge_compress_tar {
    	input:
    		fastqs1 = sra_fastq.fastq1,
    		fastqs2 = sra_fastq.fastq2,
            title = title
    }
	output{
		File tar_fastq = merge_compress_tar.tar_fastq
    }
}
task make_text {
	String buckets
    String title
    command <<<
    echo ${buckets} | tr " " "\n" > temp.txt
        >>>
    output {
    	Array[File] temp = read_lines("temp.txt")
    }
    runtime {
        disks: 			"local-disk 1 SSD"
        cpu: 			1
        memory: 		"1 GB"
        docker: 		"litd/docker-cellranger"
        preemptible: 	1
    }
}
task srr_id {
	String bucket
    
    command <<<
    SRRID=`echo ${bucket} | cut -d/ -f5`
    echo $SRRID
    
        >>>
    output {
    	String id = read_string(stdout())
    }
    runtime {
        disks: 			"local-disk 1 SSD"
        cpu: 			1
        memory: 		"1 GB"
        docker: 		"litd/docker-cellranger"
        preemptible: 	1
    }
}

 task sra_fastq {
    File bucket
    String srr
    command <<<
    
    # Convert to FASTQ
    
      fastq-dump --split-files -O . ${bucket}
		
        >>>
    output {
    	File fastq1 = "${srr}_1.fastq"
    	File fastq2 = "${srr}_2.fastq"
    }
    runtime {
        disks: 			"local-disk 400 HDD"
        cpu: 			1
        memory: 		"8 GB"
        docker: 		"flemoine/sratoolkit"
        preemptible: 	3
    }
}
 task merge_compress_tar {
    Array[File] fastqs1
    Array[File] fastqs2
    String title
    
    command <<<
    
    cat ${sep=' ' fastqs1} | pigz -p4 > ${title}"_S1_L001_R1_001.fastq.gz"
    cat ${sep=' ' fastqs2} | pigz -p4 > ${title}"_S1_L001_R2_001.fastq.gz"
    tar -cvf my_files.tar "${title}_S1_L001_R1_001.fastq.gz" "${title}_S1_L001_R2_001.fastq.gz"
    ls -lah
		>>>
	output {
    	File tar_fastq = "my_files.tar"
    }
    runtime {
        disks: 			"local-disk 400 HDD"
        cpu: 			4
        memory: 		"16 GB"
        docker: 		"rtibiocloud/pigz:v2.4_b243f9"
        preemptible: 	3
    }
}