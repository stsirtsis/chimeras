import sys
import csv

resFiles=list(sys.argv)[1:]

for f in resFiles:
    with open(f, 'rb') as csvfile1:
        with open('csv_'+f, 'wb') as csvfile2:
            lst = []
            spamreader=csv.reader(csvfile1, delimiter='\t', quotechar='|')
            spamwriter = csv.writer(csvfile2, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            i=0
            for row in spamreader:
                lst.append(row[2])
                i+=1
                if (i%100==0):
                    spamwriter.writerow(lst)
                    lst=[]
