class Config:
    def __init__ (self):
        self.result_dirs=[]
        self.genelist_paths=[]
        self.genomics_paths=[]
        with open ("config/config.txt",'r') as config_file:
            for line in config_file:
                    #Sanitize
                    line = line.strip().replace(" ","")

                    if "result_dir" == line.split('=')[0].strip():                
                        #Adds the trailing / in case it's forgotten in the config.txt
                        #to just avoid errors altogether
                        if line[-1] != '/':
                            line += '/'
                        self.result_dirs.append(line.split('=')[-1])
                    elif "genelist_path" == line.split('=')[0].strip():
                        self.genelist_paths.append(line.split('=')[-1])
                    elif "genomics_path" == line.split('=')[0].strip():
                        if line[-1] != '/':
                            line += '/'
                        self.genomics_paths.append(line.split('=')[-1])

    def getResultsDirs(self):
        return self.result_dirs

    def getGenelistPaths(self):
        return self.genelist_paths

    def getGenomicsPaths(self):
        return self.genomics_paths