#Didi Ren
#ask user to input the two gff files
import sys
if '-g1' in sys.argv:
    position1=sys.argv.index('-g1')
    inputf1=sys.argv[position1+1]

    

if '-g2' in sys.argv:
    position2=sys.argv.index('-g2')
    inputf2=sys.argv[position2+1]   



#take the gff file of the old jgi-version and obtain information

inf1=open(inputf1,'r')
jgimrna={}
jgicds={}
a=0
numericID=''     
cds=[]
for line in inf1:
    line=line.strip()
    fields=line.split('\t')
    
    if len(fields)>2:
        if fields[2]=='gene' and len(cds)>0:
            jgicds[ID]=cds
        subfields=fields[8].split(';')
        
        if fields[2]=='mRNA':
            a+=1
            for i in range(0,len(subfields)):
                if 'Parent=' in subfields[i]:
                    ID=subfields[i][7:]
                if 'JGI_V11_' in subfields[i]:
                    numericID=subfields[i][11:]
            #a dictionary of mrna line info connected to keys of jgi gene id
            if len(numericID)>0:
                jgimrna[ID]=(fields[0],fields[3],fields[4],numericID)
            else:
                jgimrna[ID]=(fields[0],fields[3],fields[4])

            
            cds=[]
        
        if fields[2]=='CDS':
            #cds will be a list of tuples with each CDS's coordinate in one gene
            cds.append((fields[0],fields[3],fields[4],fields[6],fields[7])) 
              
#a dictionary of jgi gene id to keys of all CDs coordinate
jgicds[ID]=cds           
inf1.close()


#take the gff file of the new-version and obtain information
inf2=open(inputf2,'r')
didimrna={}
didicds={}
b=0 
c=0
d=0
didiAED=''
didigenename=''
cds=[]
didiipr=[]
didiRegion=[]
didiiprdict={}
for line in inf2:
    line=line.strip()
    fields=line.split('\t')
    if len(fields)>7:
        subfields=fields[8].split(';')
        if fields[2]=='gene' and len(cds)>0:
            didicds[ID]=cds

        if fields[2]=='mRNA':
            b+=1
            for i in range(0,len(subfields)):
                if 'ID=' in subfields[i]:
                    ID=subfields[i][3:]
                if '_AED=' in subfields[i]:
                    didiAED=subfields[i][5:]
                if 'Note=' in subfields[i]:
                    didigenename=subfields[i][5:]
                if 'Dbxref' in subfields[i]:
                    d+=1
                    ssubfields=subfields[i].split(',')
                    for i in range(0,len(ssubfields)):
                        if 'InterPro' in ssubfields[i]:
                            sssubfields=ssubfields[i].split(':')
                            didiipr.append(sssubfields[1])

            #a dictionary of mrna line info connected to keys of didi gene id
            didimrna[ID]=(fields[0],fields[3],fields[4],didiAED,didigenename)
            didiiprdict[ID]=didiipr
            if len(didiipr)>0:
                c+=1
            cds=[]
            didiipr=[]
        if fields[2]=='CDS':
            #cds will be a list of tuples with each CDS's coordinate in one gene
            cds.append((fields[0],fields[3],fields[4],fields[6],fields[7]))
    #get all info for the regions in 2017version that do have RNA evidence but not predict a gene
    if len(fields)>2:
        if len(fields[1])>2:
            didiRegion.append((fields[0],fields[3],fields[4],fields[5],fields[6]))    
#a dictionary of jgi gene id to keys of all CDs coordinate
didicds[ID]=cds
inf2.close()

#because it is impossible for old gff file to be updated with interproscan resutls directly
#Here, this step take the file from old jgi-version interproscan resutls
if '-interpro' in sys.argv:
    position3=sys.argv.index('-interpro')
    inputf3=sys.argv[position3+1]  
    inf29=open(inputf3,'r')

    jgiipr={}
    iprlist=[]

    l1=inf29.readline()
    f1=l1.split('\t')
    jgiID=f1[0]
    e=2
    for line in inf29:
        line=line.strip()
        fields=line.split('\t')

        #for my case, there are multiple IDs in the first column
        if len(fields)>12:
            if 'IPR' in fields[11]:

                if fields[0] == jgiID:
                    if fields[11] not in iprlist:
                        iprlist.append(fields[11])
                else:
                    jgiipr[jgiID]=iprlist
                    jgiID=fields[0]
                    e+=1
                    iprlist=[]
                    iprlist.append(fields[11])
    jgiipr[jgiID]=iprlist
    p=[]
    inf29.close()

if '-gl' in sys.argv:
    position4=sys.argv.index('-gl')
    inputf4=sys.argv[position4+1] 

    jgiAED={}
    inf3=open(inputf4,'r')
    for line in inf3:
        line=line.strip()
        fields=line.split('\t')
        if len(fields)>8:
            subfields=fields[8].split(';')

            if fields[2]=='mRNA':
                for i in range(0,len(subfields)):
                    if 'Parent=' in subfields[i]:
                        ID=subfields[i][7:]
                    if '_AED' in subfields[i]:
                        AEDscore=subfields[i][5:]
                        jgiAED[ID]=AEDscore
    inf3.close()
    print('jgi_genenumber: ',a,'\n','jgi_genenumber.have.conserved.domain.IPR/GPID: ',e,'\n','didi_genenumber: ',b,'\n','didi_genenumber.havehits.domain.database: ',d,'\n','didi_genenumber.have.conserved.domain.IPR/GPID: ',c,'\n')




l=0
z=0
jgimatchgenes=[]
jgi1cdsgenes=[]
#obtain genes from JGI that have exact CDS with genes of didi-version
for k,v in jgicds.items():
    for k1,v1 in didicds.items():
        if set(v)==set(v1):
            l+=1
            jgimatchgenes.append((k,k1))
#obtain genes from JGI that have one least one exact CDS with genes of didi-version
               
        else:
            for x in range(0,len(v)):
                if v[x] in v1:
                    jgi1cdsgenes.append((k,k1))

jgimatchgenes_noreplicate=list(set(jgimatchgenes))               
print('JGI_genenumber.match.didiverison: ', len(jgimatchgenes), '\n')

jgi1cdsgenes_noreplicate=list(set(jgi1cdsgenes))
print('JGI_genenumber.atleast.1cds.match.didiverison: ', len(jgi1cdsgenes_noreplicate), '\n')



collect=[]
#obtain the JGI genes that share no exact CDS with didi-version
for i in jgimatchgenes_noreplicate:
    collect.append(i[0])
for j in jgi1cdsgenes_noreplicate:
    collect.append(j[0])
print('JGI_genenumber.atleast.1cds.match.multiple.hits: ', len(jgi1cdsgenes_noreplicate)-(len(list(set(collect)))-len(jgimatchgenes)), '\n')

nomatch=[]    
for k,v in jgimrna.items():
    if k not in collect:
        #nomatch is a list of tuple containing each non-match gene's mrna coordinate and name
            nomatch.append((k,v[0],v[1],v[2]))
       
print('JGI_genenumber.no.cdsmatch.didiversion: ', len(nomatch), '\n')

nomatch_overlap=[]
for i in nomatch:
    for k,v1 in didimrna.items():
        #if the scaffold info of a mrna from 2009version equals to any mrna's scaffold info from 2017version
        if i[1]==v1[0]:
            #if the locus of this mrna from 2009version is overlaping with any mrna' locus from 2017version
            if float(i[2])<float(v1[2]) and float(i[3])>float(v1[1]):
                #t is the list of the gene's name from 2009version that satisfy the above two conditions
                nomatch_overlap.append((i[0],k))
nomatch_overlap_noreplicate=list(set(nomatch_overlap))
print('JGI_genenumber.nomatch.cds.overlaps.didiversion: ', len(nomatch_overlap_noreplicate), '\n')


nomatch_nooverlap=[]
n=[]

for i in nomatch_overlap_noreplicate:
    n.append(i[0])
for j in nomatch:
    if j[0] not in n:
        #m is the list of the gene's name from 2009 version that no overlap with any gene from 2017 version
        nomatch_nooverlap.append(j[0])
print('JGI_genenumber.nomatch.cds.no.overlaps.didiversion: ', len(nomatch_nooverlap), '\n')


#output the MATCH file that share same IPR domains

q=0
outf1=open('MATCH.txt','w')
for i in jgimatchgenes_noreplicate:
    q+=1
    outf1.write(i[0]+'\t')
    for a in range(0,len(jgimrna[i[0]])):
        outf1.write(jgimrna[i[0]][a]+'\t')
    outf1.write(jgiAED[i[0]]+'\t')
    if i[0] in jgiipr.keys():
        outf1.write(str(jgiipr[i[0]])+'\t')
    else:
        outf1.write('[]'+'\t')

    outf1.write(i[1]+'\t')
    for b in range(0,len(didimrna[i[1]])):
        outf1.write(str(didimrna[i[1]][b])+'\t')
    outf1.write(str(didiiprdict[i[1]])+'\t')
    outf1.write('\n')

outf1.close()

#output the similar file that share same IPR domains
#output the different file that share some different IPR domains
j=0
p=0
k=0
f=0
r=0
y=0
t=0
outf2=open('SIMILAR.txt','w')
outf3=open('DIFFERENT.txt','w')
for i in jgi1cdsgenes_noreplicate:
    if i[0] in jgiipr.keys():
        if set(jgiipr[i[0]])==set(didiiprdict[i[1]]):
            j+=1
            outf2.write(i[0]+'\t')
            for a in range(0,len(jgimrna[i[0]])):
                outf2.write(jgimrna[i[0]][a]+'\t')
            outf2.write(jgiAED[i[0]]+'\t')
            outf2.write(str(jgiipr[i[0]])+'\t')
            outf2.write(i[1]+'\t')
            for b in range(0,len(didimrna[i[1]])):
                outf2.write(str(didimrna[i[1]][b])+'\t')
            outf2.write(str(didiiprdict[i[1]])+'\t')
            outf2.write('\n')
        else:
            p+=1
            outf3.write(i[0]+'\t')
            for a in range(0,len(jgimrna[i[0]])):
                outf3.write(jgimrna[i[0]][a]+'\t')
            outf3.write(jgiAED[i[0]]+'\t')
            outf3.write(str(jgiipr[i[0]])+'\t')
            outf3.write(i[1]+'\t')
            for b in range(0,len(didimrna[i[1]])):
                outf3.write(str(didimrna[i[1]][b])+'\t')
            outf3.write(str(didiiprdict[i[1]])+'\t')
            outf3.write('\n')
    elif len(didiiprdict[i[1]])==0:
        k+=1
        outf2.write(i[0]+'\t')
        for a in range(0,len(jgimrna[i[0]])):
            outf2.write(jgimrna[i[0]][a]+'\t')
        outf2.write(jgiAED[i[0]]+'\t')
        outf2.write('[]'+'\t')
        outf2.write(i[1]+'\t')
        for b in range(0,len(didimrna[i[1]])):
            outf2.write(str(didimrna[i[1]][b])+'\t')
        outf2.write(str(didiiprdict[i[1]])+'\t')
        outf2.write('\n')
    elif len(didiiprdict[i[1]])>0:
        f+=1
        outf3.write(i[0]+'\t')
        for a in range(0,len(jgimrna[i[0]])):
            outf3.write(jgimrna[i[0]][a]+'\t')
        outf3.write(jgiAED[i[0]]+'\t')
        outf3.write('[]'+'\t')
        outf3.write(i[1]+'\t')
        for b in range(0,len(didimrna[i[1]])):
            outf3.write(str(didimrna[i[1]][b])+'\t')
        outf3.write(str(didiiprdict[i[1]])+'\t')
        outf3.write('\n')
for i in nomatch_overlap_noreplicate:
    if i[0] in jgiipr:
        if set(jgiipr[i[0]])==set(didiiprdict[i[1]]):
            r+=1
            outf2.write(i[0]+'\t')
            for a in range(0,len(jgimrna[i[0]])):
                outf2.write(jgimrna[i[0]][a]+'\t')
            outf2.write(jgiAED[i[0]]+'\t')
            outf2.write(str(jgiipr[i[0]])+'\t')
            outf2.write(i[1]+'\t')
            for b in range(0,len(didimrna[i[1]])):
                outf2.write(str(didimrna[i[1]][b])+'\t')
            outf2.write(str(didiiprdict[i[1]])+'\t')
            outf2.write('\n')
        else:
            y+=1
            outf3.write(i[0]+'\t')
            for a in range(0,len(jgimrna[i[0]])):
                outf3.write(jgimrna[i[0]][a]+'\t')
            outf3.write(jgiAED[i[0]]+'\t')
            outf3.write(str(jgiipr[i[0]])+'\t')
            outf3.write(i[1]+'\t')
            for b in range(0,len(didimrna[i[1]])):
                outf3.write(str(didimrna[i[1]][b])+'\t')
            outf3.write(str(didiiprdict[i[1]])+'\t')
            outf3.write('\n')
    else:
        t+=1
        outf3.write(i[0]+'\t')
        for a in range(0,len(jgimrna[i[0]])):
            outf3.write(jgimrna[i[0]][a]+'\t')
        outf3.write(jgiAED[i[0]]+'\t')
        outf3.write('[]'+'\t')
        outf3.write(i[1]+'\t')
        for b in range(0,len(didimrna[i[1]])):
            outf3.write(str(didimrna[i[1]][b])+'\t')
        outf3.write(str(didiiprdict[i[1]])+'\t')
        outf3.write('\n')  
outf2.close()
outf3.close()




x=[]

for l in didiRegion:
    for i in nomatch_nooverlap:
        #if the region that predict no gene but with RNA-reads in 2017version also overlap with any gene in 2009version
        if jgimrna[i][0]==l[0]:
            if float(jgimrna[i][1])<float(l[2]) and float(jgimrna[i][2])>float(l[1]):
                x.append(i)
                
print('JGI_genenumber.nomatch.cds.no.overlaps.with.mrna.support.didiversion: ', len(list(set(x))), '\n')           

#for any genes in 2009version do not match or overlap with any gene in 2017version, nor overlap with any RNA-reads-no-gene region in 2017version
#indicating no RNA were transcribed from this predicted gene in 2009version
g=0
h=0
outf=open('NO_GENE.txt','w')
outf4=open('NO_GENE_WITH_RNA.txt','w')

for i in nomatch_nooverlap:
    if i not in list(set(x)):
        outf.write(i+'\t')
        for a in range(0,len(jgimrna[i])):
            outf.write(jgimrna[i][a]+'\t')
        outf.write(jgiAED[i]+'\t')
        outf.write('\n')
        g+=1
    else:
        h+=1
        outf4.write(i+'\t')
        for a in range(0,len(jgimrna[i])):
            outf4.write(jgimrna[i][a]+'\t')
        outf4.write(jgiAED[i]+'\t')
        outf4.write('\n')

outf.close()
outf4.close()
print('NO_GENE.txt genenumber: ', g, '\n')
print('NO_GENE.RNA.support.txt genenumber: ', h, '\n')           

print('SIMILAR.txt genenumber: ',j+k+r,'\n','   In them, atleast.1cds.match.domain.match: ',j,'\n','    In them, atleast.1cds.match.both.nodomains:',k,'\n', '  In them, nomatch.cds.overlap.domain.match:',r,'\n')           
print('DIFFERENT.txt genenumber: ',p+f+y+t,'\n','   In them, atleast.1cds.match.domain.not-match: ',p,'\n','    In them, atleast.1cds.match.jgi.nodomain.didi.withdomains:',f, '\n','   In them, nomatch.cds.overlap.domain.not-match:',y,'\n', '   In them, nomatch.cds.overlap.jgi.withoutdomains:',t,'\n')           
print('MATCH.txt genenumber: ',q, '\n')           










       

