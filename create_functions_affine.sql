create or replace function needleman_affine(ids integer, sequences text, score_match integer, score_mismatch integer, score_gap_open integer, score_gap_extend integer) 
returns table(id1 integer, id2 integer, align1 text, align2 text, score integer) language python {
    from itertools import product
    import itertools
    from collections import deque
    import math
    
    id1, id2, align1, align2, score = [], [], [], [], []

    score_match, score_mismatch, score_gap_open, score_gap_extend = score_match[0], score_mismatch[0], score_gap_open[0], score_gap_extend[0]
    
    for (p,x), (q,y) in list(itertools.combinations(enumerate(sequences), 2)):
      N, M = len(x), len(y)
      
      s = lambda a, b: score_match if a==b else score_mismatch

      DIAG = -1, -1
      LEFT = -1, 0
      UP = 0, -1
      
      # Create tables F, X, Y, V and Ptr
      F, X, Y, V = {}, {}, {}, {}
      Ptr = {}
      
      F[-1, -1], X[-1,-1], Y[-1,-1], V[-1,-1] = 0, score_gap_open + score_gap_extend, score_gap_open + score_gap_extend, 0
      for i in range(N):
        F[i, -1] = -math.inf
        X[i, -1] = score_gap_open + (i+1) * score_gap_extend
        Y[i, -1] = -math.inf
        V[i, -1] = score_gap_open + (i+1) * score_gap_extend
      for j in range(M):
        F[-1, j] = -math.inf
        X[-1, j] = -math.inf
        Y[-1, j] = score_gap_open + (j+1) * score_gap_extend
        V[-1, j] = score_gap_open + (j+1) * score_gap_extend
        
      
      option_Ptr = DIAG, LEFT, UP
      for i, j in product(range(N), range(M)):
        F[i, j] = V[i - 1, j - 1] + s(x[i], y[j])
        
        X[i,j] = max(F[i-1,j] + score_gap_open + score_gap_extend, X[i-1,j] + score_gap_extend)
        Y[i,j] = max(F[i,j-1] + score_gap_open + score_gap_extend, Y[i,j-1] + score_gap_extend)
        
        option_V = (F[i, j], X[i, j], Y[i, j])
        V[i, j], Ptr[i, j] = max(zip(option_V, option_Ptr), key = lambda k: k[0])
       
      # Work backwards from (N - 1, M - 1) to (0, 0)
      # to find the best alignment.
      alignment = deque()
      i, j = N - 1, M - 1
      while i >= 0 and j >= 0:
        direction = Ptr[i, j]
        if direction == DIAG:
            element = i, j
        elif direction == LEFT:
            element = i, None
        elif direction == UP:
            element = None, j
        alignment.appendleft(element)
        di, dj = direction
        i, j = i + di, j + dj
      while i >= 0:
        alignment.appendleft((i, None))
        i -= 1
      while j >= 0:
        alignment.appendleft((None, j))
        j -= 1
      
      first = "".join(
        "-" if i is None else x[i] for i, _ in alignment
      )
      second = "".join(
        "-" if j is None else y[j] for _, j in alignment
      )
      
      curr_score = V[N-1, M-1]
              
      id1.append(ids[p])
      id2.append(ids[q])
      align1.append(first)
      align2.append(second)
      score.append(curr_score)
      
    
    result = dict()
    result['id1'] = id1
    result['id2'] = id2
    result['align1'] = align1
    result['align2'] = align2
    result['score'] = score
    
    return result
};

create or replace function text_align_affine(seq1 text, seq2 text, score_match integer, score_mismatch integer, score_gap_open integer, score_gap_extend integer) 
returns table(align1 text, align2 text, score integer) language python {
    from itertools import product
    import itertools
    from collections import deque
    import math
    
    align1, align2, score = [], [], []

    x, y = seq1, seq2
    N, M = len(x), len(y)

    score_match, score_mismatch, score_gap_open, score_gap_extend = score_match, score_mismatch, score_gap_open, score_gap_extend
      
    s = lambda a, b: score_match if a==b else score_mismatch

    DIAG = -1, -1
    LEFT = -1, 0
    UP = 0, -1
    
    # Create tables F, X, Y, V and Ptr
    F, X, Y, V = {}, {}, {}, {}
    Ptr = {}
    
    F[-1, -1], X[-1,-1], Y[-1,-1], V[-1,-1] = 0, score_gap_open + score_gap_extend, score_gap_open + score_gap_extend, 0
    for i in range(N):
      F[i, -1] = -math.inf
      X[i, -1] = score_gap_open + (i+1) * score_gap_extend
      Y[i, -1] = -math.inf
      V[i, -1] = score_gap_open + (i+1) * score_gap_extend
    for j in range(M):
      F[-1, j] = -math.inf
      X[-1, j] = -math.inf
      Y[-1, j] = score_gap_open + (j+1) * score_gap_extend
      V[-1, j] = score_gap_open + (j+1) * score_gap_extend
      
    
    option_Ptr = DIAG, LEFT, UP
    for i, j in product(range(N), range(M)):
      F[i, j] = V[i - 1, j - 1] + s(x[i], y[j])
      
      X[i,j] = max(F[i-1,j] + score_gap_open + score_gap_extend, X[i-1,j] + score_gap_extend)
      Y[i,j] = max(F[i,j-1] + score_gap_open + score_gap_extend, Y[i,j-1] + score_gap_extend)
      
      option_V = (F[i, j], X[i, j], Y[i, j])
      V[i, j], Ptr[i, j] = max(zip(option_V, option_Ptr), key = lambda k: k[0])
     
    # Work backwards from (N - 1, M - 1) to (0, 0)
    # to find the best alignment.
    alignment = deque()
    i, j = N - 1, M - 1
    while i >= 0 and j >= 0:
      direction = Ptr[i, j]
      if direction == DIAG:
          element = i, j
      elif direction == LEFT:
          element = i, None
      elif direction == UP:
          element = None, j
      alignment.appendleft(element)
      di, dj = direction
      i, j = i + di, j + dj
    while i >= 0:
      alignment.appendleft((i, None))
      i -= 1
    while j >= 0:
      alignment.appendleft((None, j))
      j -= 1
    
    first = "".join(
      "-" if i is None else x[i] for i, _ in alignment
    )
    second = "".join(
      "-" if j is None else y[j] for _, j in alignment
    )
    
    curr_score = V[N-1, M-1]
            
    align1.append(first)
    align2.append(second)
    score.append(curr_score)
      
    
    result = dict()
    result['align1'] = align1
    result['align2'] = align2
    result['score'] = score
    
    return result
};

create or replace function create_lookup_affine(align1 text, align2 text, score_match integer, score_mismatch integer, score_gap_open integer, score_gap_extend integer, k integer) 
returns table(align1 text, align2 text, score integer) language python {
  import math
  
  k = k[0]
  threshold = math.ceil(k/2)
  kmers1, kmers2, score = [], [], [] 
  
  for i in range(len(align1)):
    al1, al2 = align1[i], align2[i]
    
    for j in range(len(al1)-k+1):
      kmer1 = al1[j:j+k]
      kmer2 = al2[j:j+k]
      
      curr_score = 0
      first_gap = True
      for p in range(k):
        if al1[j+p] == al2[j+p]:
          curr_score += score_match[0]
          first_gap = True
        elif al1[j+p] == '-' or al2[j+p] == '-':
          if first_gap:
            curr_score += score_gap_open[0]
            first_gap = False
          else:
            curr_score += score_gap_extend[0]
        else:
          curr_score += score_mismatch[0]
          first_gap = True
        
      if curr_score >= threshold:
        kmers1.append(kmer1)
        kmers2.append(kmer2)
        score.append(curr_score)
        
  return [kmers1, kmers2, score]

};

create or replace procedure update_lookup_affine(score_match integer, score_mismatch integer, score_gap_open integer, score_gap_extend integer, k integer, samplesize integer, upperlimit integer)
begin
  truncate table aligned_subseq CASCADE;
  truncate table lookup CASCADE;
  truncate table sample_seq;
  
  insert into sample_seq (SELECT * FROM seq_sample(samplesize, upperlimit));
  
  insert into lookup(align1, align2, score)
  select * from create_lookup((select align1, align2, score_match, score_mismatch, score_gap_open,score_gap_extend, k 
  from needleman_affine((select id, seq, score_match, score_mismatch, score_gap_open, score_gap_extend
  from sample_seq ))));

  insert into aligned_subseq(seq) 
  select * from (
  (select replace(align1, '-', '') as "s" from lookup)
  union distinct 
  (select replace(align2, '-', '') as "s" from lookup)
  ) as u
  where not exists (select x.seq from aligned_subseq x where x.seq=u.s);
  
  update lookup set seqid1 = (select id from aligned_subseq where aligned_subseq.seq = replace(align1, '-', '')),
  seqid2 = (select id from aligned_subseq where aligned_subseq.seq = replace(align2, '-', ''));
end;

create or replace function use_subalign_affine(ids integer, align1 text, align2 text, score integer, ind1 integer, ind2 integer, first integer, s1 text, s2 text, score_match integer, score_mismatch integer, score_gap_open integer, score_gap_extend integer)
returns text
language python {
from itertools import product
import itertools
from collections import deque
import datetime
import math

isLogging = False

if isLogging:
  log = open(r"monet_affine.log", 'a')
  log.write("\nSTARTING ALIGNMENT.. " + str(datetime.datetime.now()) + '\n')

def align(x,y):  
    N, M = len(x), len(y)
    
    s = lambda a, b: score_match if a==b else score_mismatch

    DIAG = -1, -1
    LEFT = -1, 0
    UP = 0, -1

    # Create tables F and Ptr
    F, X, Y, V = {}, {}, {}, {}
    Ptr = {}

    F[-1, -1], X[-1,-1], Y[-1,-1], V[-1,-1] = 0, score_gap_open, score_gap_open, 0
    for i in range(N):
      F[i, -1] = -math.inf
      X[i, -1] = score_gap_open + (i+1) * score_gap_extend
      Y[i, -1] = -math.inf
      V[i, -1] = score_gap_open + (i+1) * score_gap_extend
    for j in range(M):
      F[-1, j] = -math.inf
      X[-1, j] = -math.inf
      Y[-1, j] = score_gap_open + (j+1) * score_gap_extend
      V[-1, j] = score_gap_open + (j+1) * score_gap_extend
      
    
    option_Ptr = DIAG, LEFT, UP
    for i, j in product(range(N), range(M)):
      F[i, j] = V[i - 1, j - 1] + s(x[i], y[j])
      
      X[i,j] = max(F[i-1,j] + score_gap_open, X[i-1,j] + score_gap_extend)
      Y[i,j] = max(F[i,j-1] + score_gap_open, Y[i,j-1] + score_gap_extend)
      
      option_V = (F[i, j], X[i, j], Y[i, j])
      V[i, j], Ptr[i, j] = max(zip(option_V, option_Ptr), key = lambda k: k[0])
     
    # Work backwards from (N - 1, M - 1) to (0, 0)
    # to find the best alignment.
    alignment = deque()
    i, j = N - 1, M - 1
    while i >= 0 and j >= 0:
      direction = Ptr[i, j]
      if direction == DIAG:
          element = i, j
      elif direction == LEFT:
          element = i, None
      elif direction == UP:
          element = None, j
      alignment.appendleft(element)
      di, dj = direction
      i, j = i + di, j + dj
    while i >= 0:
      alignment.appendleft((i, None))
      i -= 1
    while j >= 0:
      alignment.appendleft((None, j))
      j -= 1
    
    first = "".join(
      "-" if i is None else x[i] for i, _ in alignment
    )
    second = "".join(
      "-" if j is None else y[j] for _, j in alignment
    )

    curr_score = V[N-1, M-1]
  
    if isLogging:
      log.write('SUBALIGNMENT: \n ' + first + '\n ' + second + '\n')
    
    return first, second, curr_score

def check_intersect(arr, ind, r):
    summa = sum(arr[ind:ind+r])
    if summa > 0:
      return False
    else:
      arr[ind:ind+r] = [1] * r
      return True
try:     
  seq1 = list(s1)
  seq2 = list(s2)
  
  used1 = [0] * len(seq1)
  used2 = [0] * len(seq2)


  for i in range(len(ids)):
      k = len(align1[i])
      
      if first[i] == 1:
        curr_ind1, curr_ind2 = ind1[i]-1, ind2[i]-1 #IND WITH NO ZERO
        curr_align1, curr_align2 = align1[i], align2[i]
      else:
        curr_ind1, curr_ind2 = ind1[i]-1, ind2[i]-1 #IND WITH NO ZERO
        curr_align1, curr_align2 = align2[i], align1[i]
      
      curr_al_len1, curr_al_len2 = len(curr_align1.replace('-', '')), len(curr_align2.replace('-', ''))

      if check_intersect(used1, curr_ind1, curr_al_len1) and check_intersect(used2, curr_ind2, curr_al_len2):
        if isLogging:
          log.write('INSERTING: ' + curr_align1 + ', ' + curr_align2 + ' at ' + str(curr_ind1) + ', ' + str(curr_ind2) + '\n')
  
        seq1[curr_ind1] = '[' + curr_align1 + "#" + str(score[i]) + ']'
        seq2[curr_ind2] = '[' + curr_align2 + ']'
        
        seq1[curr_ind1+1:curr_ind1+curr_al_len1] = ['_'] * (curr_al_len1-1)
        seq2[curr_ind2+1:curr_ind2+curr_al_len2] = ['_'] * (curr_al_len2-1)
            
        if isLogging:
          log.write('RESULT: \n ' + ''.join(seq1) + '\n ' + ''.join(seq2) + '\n')
  
  #remove used bases:
  seq1 = ''.join(seq1).replace('_','')
  seq2 = ''.join(seq2).replace('_','')
  
  if isLogging:
    log.write('ALL SUBALIGNMENTS INSERTED: \n ' + seq1 + '\n ' + seq2 + '\n')
  
  res1, res2, score = '', '', 0
  
  while '[' in seq1:
      i1 = seq1.find('[') 
      i2 = seq2.find('[')
      j1 = seq1.find(']') 
      j2 = seq2.find(']')
      
      if len(seq1[:i1]) > 0 or len(seq2[:i2]) > 0:
        first, sec, currscore = align(seq1[:i1], seq2[:i2])
      else:
        first, sec, currscore = '', '', 0
      
      tmp = seq1[i1+1:j1].split("#")
      
      res1 += first + tmp[0]
      res2 += sec + seq2[i2+1:j2]
      score += currscore + int(tmp[1])
      
      seq1, seq2 = seq1[j1+1:], seq2[j2+1:]
      
      if isLogging and (len(seq1) > 0 or len(seq2) > 0):
        log.write('REMAINING : \n ' + seq1 + ' \n ' + seq2 + '\n')
  
  if len(seq1) > 0 or len(seq2) > 0:
    first, sec, currscore = align(seq1, seq2)
    res1 += first
    res2 += sec
    score += currscore

  if isLogging:
    log.write('RESULT : \n ' + res1 + '  \n ' + res2 + ' \n SCORE: ' + str(score) + '\n')
    log.write('==================================================')
    log.close()
  
  return str(score) + " , " + res1 +" , " + res2
  
except Exception as e:
  if isLogging:
    log.write("EXCEPTION\n")
    log.write(str(e))
    log.close()
  
  return 'ERROR see log'
};

create or replace function use_lookup_affine(id1 integer, id2 integer, score_match integer, score_mismatch integer, score_gap integer)
returns text
begin
  declare result text, seq1 text, seq2 text;
  
  select seq into seq1 from seq_data where id = id1;
  select seq into seq2 from seq_data where id = id2;
  
  select use_subalign_affine(id,align1,align2,score,ind1,ind2,first, seq1, seq2, score_match, score_mismatch, score_gap_open, score_gap_extend) into result
  from get_common_subaligns(id1,id2) where -(abs(ind1-ind2)*score_gap_extend + score_gap_open) < score;
    
  if (result is null)
    THEN select concat(concat(concat(concat(align1,','),align2),','),score) into result from 
      text_align_affine(seq1, seq2, score_match, score_mismatch, score_gap_open, score_gap_extend);
  end if;
  
  return result;
end;

create or replace function use_lookup_simple_affine(id1 integer, id2 integer, score_match integer, score_mismatch integer, score_gap_open integer, score_gap_extend integer)
returns text
begin
  declare result text, seq1 text, seq2 text;
  
  select seq into seq1 from seq_data where id = id1;
  select seq into seq2 from seq_data where id = id2;
  
  select use_subalign_affine(id,align1,align2,score,ind1,ind2,first, seq1, seq2, score_match, score_mismatch, score_gap_open, score_gap_extend) into result
  from get_common_subaligns(id1,id2) where -(abs(ind1-ind2)*score_gap_extend + score_gap_open) < score;
  
  return result;
end;