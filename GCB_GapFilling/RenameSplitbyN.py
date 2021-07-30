import sys

faFILE = sys.argv[1]

scaf = ''
out = faFILE.split(".fa")[0] + "_SplitbyN.fa"
outFile = open(out, "w")


with open(faFILE, "r") as file:
	for line in file:
		if ">" in line:
			if scaf != '':
				c = 0
                                scaf_list = scaf.split("N"*20)
                 		for i in range(len(scaf_list)):
			                if len(scaf_list[i]) > 0:
						if 'N' in scaf_list[i][:19]:
							t = scaf_list[i][:19].replace("N", '')
							scaf_list[i] = t + scaf_list[i][19:]
                			        scafLen = len(scaf_list[i])
						outFile.write(scaf_name+"_p"+str(c)+"_length_"+str(scafLen)+"\n")
			                        cnt = 0
                		        	while cnt < scafLen:
                                			end = scafLen if cnt + 100 > scafLen else cnt + 100
		                        	        outFile.write(str(scaf_list[i][cnt:end]) + "\n")
                			                cnt += 100
        	                		c += 1

			scaf_name = line.split()[0]
			scaf = ''
			continue

		scaf += line.strip().upper()
		

if scaf != '':
	c = 0
	scaf_list = scaf.split("N"*20)
	for i in range(len(scaf_list)):
		if len(scaf_list[i]) > 0:
			if 'N' in scaf_list[i][:19]:
				t = scaf_list[i][:19].replace("N", '')
                                scaf_list[i] = t + scaf_list[i][19:]
			scafLen = len(scaf_list[i])
			outFile.write(scaf_name+"_p"+str(c)+"_length_"+str(scafLen)+"\n")
			cnt = 0
        		while cnt < scafLen:
				end = scafLen if cnt + 100 > scafLen else cnt + 100
				outFile.write(str(scaf_list[i][cnt:end]) + "\n")
				cnt += 100
			c += 1
	
#	scafLen = len(scaf)
	#outFile.write(scaf_name+"_length_"+str(scafLen)+"\n")
	#cnt = 0
	#while cnt < scafLen:
#		end = scafLen if cnt + 100 > scafLen else cnt + 100
#		outFile.write(str(scaf[cnt:end]) + "\n")
#		cnt += 100
