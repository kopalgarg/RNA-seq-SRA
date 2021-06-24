workflow SRA_Downloads{
	File genomefasta
    File genomeGTF
    Int readLength
    File genomefastaIndex

	call create_index {
    	input:
        	genomefasta = genomefasta,
            genomeGTF = genomeGTF,
            readLength = readLength,
            genomefastaIndex = genomefastaIndex
    }

}

task create_index {

	File genomefasta
    File genomeGTF
    Int readLength
    File genomefastaIndex
    
    command <<<
        
    # create index (prefix: 'index')

        mkdir out && cd out
        rsem-prepare-reference --gtf ${genomeGTF} -p 16 --star  --gtf ${genomeGTF} ${genomefasta} index --star-sjdboverhang ${readLength}
	
    	STAR --runMode genomeGenerate \
    	 	--runThreadN 16  \
         	--genomeFastaFiles ${genomefasta} \
         	--sjdbGTFfile ${genomeGTF} \
         	--sjdbOverhang ${readLength} \
         	--outFileNamePrefix index
        
        la -lah
        
        cd .. && tar -cvzf out.tar.gz out

        >>>
        
        
    output {
    
    	File rsem_index = "out.tar.gz"
    }
    
    runtime {
        disks: 			"local-disk 125 HDD"
        cpu: 			16
        memory: 		"60 GB"
        docker: 		"gargk/sra"
        preemptible: 	3
    }
}