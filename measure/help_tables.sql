CREATE TABLE alignments_help(
  id1 integer, 
  id2 integer,
  score integer,
  time double,

  constraint ah_pk primary key (id1,id2)
);

CREATE TABLE alignments_help_affine(
  id1 integer, 
  id2 integer,
  score integer,
  time double,

  constraint ahf_pk primary key (id1,id2)
);