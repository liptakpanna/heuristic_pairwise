CREATE TABLE alignments_help(
  id1 integer, 
  id2 integer,
  score integer,
  time integer
);

create or replace function needleman_affine_help(i1 integer, i2 integer, sequences text, score_match integer, score_mismatch integer, score_gap_open integer, score_gap_extend integer) 
returns text
begin
    declare s integer, time integer, start timestamp, end timestamp;

    select score into s, time into time from alignments_help where (id1 = i1 and id2 = i2) or (id2 = i1 and id1 = i2);

    if score is null then
        select now() into start;

        select score into s from needleman_affine(select id, seq, score_match, score_mismatch,score_gap_open,score_gap_extend,
        from seq_data where id in (id1,  id2));

        select now() into end;
    end if;

    insert into alignments_help values (id1, id2, s, end-start);

    return concatenate(concatenate(s, "#"), end-now);
end;
