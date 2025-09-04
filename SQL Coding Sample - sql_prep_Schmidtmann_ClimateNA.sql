/* Convert 2001 Schmidtmann data into ClimateNA input format */
CREATE TABLE "Schmidtmann_data_final_only_2001_ClimateNA" AS
SELECT
id AS ID1,
farmID AS ID2,
latitude AS lat,
longitude AS long
FROM Schmidtmann_data_final_only_2001
;

ALTER TABLE Schmidtmann_data_final_only_2001_ClimateNA
ADD COLUMN el REAL
;

UPDATE Schmidtmann_data_final_only_2001_ClimateNA
SET el = "."
;

/* Convert 2002 Schmidtmann data into ClimateNA input format */
CREATE TABLE "Schmidtmann_data_final_only_2002_ClimateNA" AS
SELECT
id AS ID1,
farmID AS ID2,
latitude AS lat,
longitude AS long
FROM Schmidtmann_data_final_only_2002
;

ALTER TABLE Schmidtmann_data_final_only_2002_ClimateNA
ADD COLUMN el REAL
;

UPDATE Schmidtmann_data_final_only_2002_ClimateNA
SET el = "."
;

/* Convert 2001 Schmidtmann data (including previously NA year values that 
were then estimated as 2001) into ClimateNA input format */
CREATE TABLE "Schmidtmann_data_final_only_2001_incl_NAs_ClimateNA" AS
SELECT
id AS ID1,
farmID AS ID2,
latitude AS lat,
longitude AS long
FROM Schmidtmann_data_final_only_2001_incl_NAs
;

ALTER TABLE Schmidtmann_data_final_only_2001_incl_NAs_ClimateNA
ADD COLUMN el REAL
;

UPDATE Schmidtmann_data_final_only_2001_incl_NAs_ClimateNA
SET el = "."
;