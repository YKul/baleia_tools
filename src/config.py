#Represents a configuration as read from a config file
#Todo: config path should actually be something that can be edited/changed
import os
class Config:
    def __init__ (self):
        #Filepaths are stored in lists. Lists are parallelized by
        #data groups defined in the configuration file.
        #Eg. result[0], genomes[0], etc all correspond to the first 
        #group defined in the config file
        self.results=[]
        self.genelists=[]
        self.genomes=[]
        self.isoforms=[]
        self.extracted=[]
        self.filtered=[]
        self.name_changed_genes=[]
        self.multispecies_genes=[]
        self.proteins=[]
        self.multispecies_proteins=[]
        self.nuc_alignments=[]
        self.cluster_jobs=[]
        with open ("config/config.txt",'r') as config_file:
            for line in config_file:
                    
                    #Sanitize
                    line = line.strip().replace(" ","")

                    if "results_dir" == line.split('=')[0].strip():                
                        #Adds the trailing / in case it's forgotten in the config.txt
                        #to just avoid errors altogether
                        if line[-1] != '/':
                            line += '/'
                        self.results.append(line.split('=')[-1])
                        os.makedirs(self.results[-1], exist_ok=True)
                    #Note that genelist is a filepath not a directory
                    elif "genelist_path" == line.split('=')[0].strip():
                        self.genelists.append(line.split('=')[-1])
                    elif "genomics_path" == line.split('=')[0].strip():
                        if line[-1] != '/':
                            line += '/'
                        self.genomes.append(line.split('=')[-1])
                        os.makedirs(self.genomes[-1], exist_ok=True)
                    elif "cluster_jobs" == line.split('=')[0].strip():
                        if line[-1] != '/':
                            line += '/'
                        self.cluster_jobs.append(line.split('=')[-1])
                        os.makedirs(self.cluster_jobs[-1], exist_ok=True)

        for result_dir in self.results:
            os.makedirs(result_dir, exist_ok=True)

            self.isoforms.append(f"{result_dir}largest-isoforms")
            os.makedirs(self.isoforms[-1], exist_ok=True)

            self.extracted.append(f"{result_dir}extracted-genes")
            os.makedirs(self.extracted[-1], exist_ok=True)

            self.filtered.append(f"{result_dir}filtered-lowq")
            os.makedirs(self.filtered[-1], exist_ok=True)

            self.name_changed_genes.append(f"{result_dir}genes-name-changed")
            os.makedirs(self.name_changed_genes[-1], exist_ok=True)

            self.multispecies_genes.append(f"{result_dir}genes-multispecies")
            os.makedirs(self.multispecies_genes[-1], exist_ok=True)

            self.proteins.append(f"{result_dir}proteins")
            os.makedirs(self.proteins[-1], exist_ok=True)

            self.multispecies_proteins.append(f"{result_dir}proteins-multispecies")
            os.makedirs(self.multispecies_proteins[-1], exist_ok=True)

            self.nuc_alignments.append(f"{result_dir}nuc_alignments")
            os.makedirs(self.nuc_alignments[-1], exist_ok=True)

        #Handles when there is only 1 shared path in the config file
        if len(self.cluster_jobs)==1 and len(self.results)>1:
            for i in range(1,len(self.results)):
                self.cluster_jobs.append(self.cluster_jobs[0])
        if len(self.genomes)==1 and len(self.results)>1:
            for i in range(1,len(self.results)):
                self.genomes.append(self.genomes[0])

    def getClusterJobsPath(self,i=None):
        if not i==None:
            return self.cluster_jobs[i]
        else:
            return self.cluster_jobs

    def getResultsPath(self,i=None):
        if not i==None:
            return self.results[i]
        else:
            return self.results

    def getGenelistsPath(self,i=None):
        if not i==None:
            return self.genelists[i]
        else:
            return self.genelists

    def getGenomesPath(self,i=None):
        if not i==None:
            return self.genomes[i]
        else:
            return self.genomes

    def getIsoformsPath(self,i=None):
        if not i==None:
            return self.isoforms[i]
        else:
            return self.isoforms

    def getExtractedPath(self,i=None):
        if not i==None:
            return self.extracted[i]
        else:
            return self.extracted

    def getFilteredPath(self,i=None):
        if not i==None:
            return self.filtered[i]
        else:
            return self.filtered

    def getNameChangedGenesPath(self,i=None):
        if not i==None:
            return self.name_changed_genes[i]
        else:
            return self.name_changed_genes

    def getMultispeciesGenesPath(self,i=None):
        if not i==None:
            return self.multispecies_genes[i]
        else:
            return self.multispecies_genes

    def getProteinsPath(self,i=None):
        if not i==None:
            return self.proteins[i]
        else:
            return self.proteins

    def getMultispeciesProteinsPath(self,i=None):
        if not i==None:
            return self.multispecies_proteins[i]
        else:
            return self.multispecies_proteins

    def getNucAlignmentsPath(self,i=None):
        if not i==None:
            return self.nuc_alignments[i]
        else:
            return self.nuc_alignments