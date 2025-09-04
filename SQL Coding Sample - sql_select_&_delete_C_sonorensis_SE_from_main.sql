/* check data by selecting rows from main occurrence table where id matches that from rows in C sonorensis southeastern US table*/
SELECT *
FROM occurrence AS main
LEFT JOIN c_sonorensis_se AS subset ON main.id = subset.id
WHERE main.id = subset.id;

/* delete rows from main occurrence table where id matches that from rows in C sonorensis southeastern US table*/
DELETE FROM occurrence AS main
WHERE EXISTS 
(SELECT * FROM c_sonorensis_se AS subset
WHERE subset.id = main.id)