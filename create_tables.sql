create table seq_data (
  id serial, 
  seq text NOT NULL
);

CREATE TABLE aligned_subseq(
  id serial, 
  seq text
);

CREATE TABLE lookup(
    id serial, 
    seqid1 int, 
    align1 text, 
    seqid2 int, 
    align2 text, 
    score integer,
    CONSTRAINT lookup_fk1 FOREIGN KEY (seqid1) REFERENCES aligned_subseq (id),
    CONSTRAINT lookup_fk2 FOREIGN KEY (seqid2) REFERENCES aligned_subseq (id)
);

create table sample_seq(
    id integer PRIMARY KEY, 
    seq text NOT NULL,
    CONSTRAINT sample_seq_fk1 FOREIGN KEY (id) REFERENCES seq_data (id)
);