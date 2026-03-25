
import os
#Remove case sensitivity
class Config:
    def __init__ (self):
        self.email=""
        self.macse=""
        self.java=""
        self.error_output=[]
        self.stream_output=[]
        self.results=[]

        with open("config/awarnach_config.txt",'r') as config_file:
            for line in config_file:                    
                    #Sanitize
                    line = line.strip().replace(" ","")
                    if "email" == line.split('=')[0].strip():                
                        #Adds the trailing / in case it's forgotten in the config.txt
                        #to just avoid errors altogether
                        if line[-1] != '/':
                            line += '/'
                        self.email=line.split('=')[-1]
                    elif "MACSE_path" == line.split('=')[0].strip():                
                        if line[-1] == '/':
                            line = line[:-1]
                        self.macse=line.split('=')[-1]
                    elif "JAVA_path" == line.split('=')[0].strip():                
                        if line[-1] == '/':
                            line = line[:-1]
                        self.java=line.split('=')[-1]
                    elif "error_output" == line.split('=')[0].strip():
                        if line[-1] != '/':
                            line += '/'
                        self.error_output.append(line.split('=')[-1])
                    elif "stream_output" == line.split('=')[0].strip():
                        if line[-1] != '/':
                            line += '/'
                        self.stream_output.append(line.split('=')[-1])
                    elif "results_dir" == line.split('=')[0].strip():                
                        #Adds the trailing / in case it's forgotten in the config.txt
                        #to just avoid errors altogether
                        if line[-1] != '/':
                            line += '/'
                        self.results.append(line.split('=')[-1])

        #Handles when there is an optional shared path in the config file
        if len(self.error_output)==1 and len(self.results)>1:
            for i in range(1,len(self.results)):
                self.error_output.append(self.error_output[0])
        if len(self.stream_output)==1 and len(self.results)>1:
            for i in range(1,len(self.results)):
                self.stream_output.append(self.stream_output[0])

    def getEmail(self):
        return self.email

    def getMACSE(self):
        return self.macse

    def getJAVA(self):
        return self.java

    def getErrorOutput(self,i=None):
        if not i==None:
            return self.error_output[i]
        else:
            return self.error_output

    def getStreamOutput(self,i=None):
        if not i==None:
            return self.stream_output[i]
        else:
            return self.stream_output

    def getResultsDirs(self,i=None):
        if not i==None:
            return self.results[i]
        else:
            return self.results