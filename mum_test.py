from pymummer import coords_file, alignment, nucmer
import pymonetdb
import time

connection = pymonetdb.connect(username="monetdb", password="monetdb", hostname="localhost", database="thesis")
cursor = connection.cursor()

start_time = time.time()

cursor.execute('SELECT * FROM seq_data where id = 1')
seq1 = cursor.fetchone()
seq1 = seq1[1]

f = open("s1.txt", "w")
f.write("> s1\n")
f.write(seq1+ '\n')
f.close()

cursor.execute('SELECT * FROM seq_data where id = 2')
seq2 = cursor.fetchone()
seq2 = seq2[1]


f = open("s2.txt", "w")
f.write("> s2\n")
f.write(seq2 + '\n')
f.close()

runner = nucmer.Runner('s1.txt', 's2.txt', 'out.txt') 
runner.run()
file_reader = coords_file.reader('out.txt')
alignments = [coord for coord in file_reader if not coord.is_self_hit()] #Remove self hits
print(alignments)

print("--- %s seconds ---" % (time.time() - start_time))