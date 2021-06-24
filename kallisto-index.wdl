workflow kallisto_index{

    File refTranscriptome

	call kallisto{
		input:
        	refTranscriptome = refTranscriptome
            
}

}

task kallisto {

	File refTranscriptome
            	
    
    command <<<
            
    	mkdir out && cd out
        kallisto index -i index ${refTranscriptome}
        la -lah
        cd .. && tar -cvzf out.tar.gz out
        

        >>>
    output {
    
        	File index = "out.tar.gz"

    	
    }

    runtime {
        disks: 			"local-disk 100 HDD"
        cpu: 			8
        memory: 		"32 GB"
        docker: 		"zlskidmore/kallisto"
        preemptible: 	3
    }
    
 }