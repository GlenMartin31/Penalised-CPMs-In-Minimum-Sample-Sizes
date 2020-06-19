SET search_path TO mimiciii;

DROP MATERIALIZED VIEW IF EXISTS Glen_Penalised_CPMs_Cohort CASCADE;

CREATE MATERIALIZED VIEW Glen_Penalised_CPMs_Cohort AS
WITH cohort AS (
SELECT ie.subject_id
, ie.hadm_id
, ie.icustay_id
, ie.intime
, ie.outtime
, ie.los AS icu_los_days
, DENSE_RANK() OVER (PARTITION BY ie.hadm_id ORDER BY ie.intime) AS icustay_seq --calculate a count of how many ICU admissions each patient had within a given hospitalisation

-- Patient level factors
, pat.dob
, ROUND((cast(ie.intime as date) - cast(pat.dob as date))/365.242, 2) AS age
, pat.gender

-- Hospital-admission-level factors
, adm.admittime
, adm.dischtime
, adm.deathtime -- provides the time of in-hospital death for the patient...only present if the patient died in-hospital
, adm.hospital_expire_flag -- indicates whether the patient died within the given hospitalization.
, DENSE_RANK() OVER (PARTITION BY adm.subject_id ORDER BY adm.admittime) AS hospstay_seq --define the number of hostpial admissions for each patient
, adm.admission_type
, CASE WHEN adm.ethnicity IN
  (
       'WHITE'
     , 'WHITE - RUSSIAN' 
     , 'WHITE - OTHER EUROPEAN' 
     , 'WHITE - BRAZILIAN'
     , 'WHITE - EASTERN EUROPEAN' 
  ) THEN 'white' --define a white ethnic group
  WHEN adm.ethnicity IN
  (
      'BLACK/AFRICAN AMERICAN' 
    , 'BLACK/CAPE VERDEAN' 
    , 'BLACK/HAITIAN' 
    , 'BLACK/AFRICAN' 
    , 'CARIBBEAN ISLAND' 
  ) THEN 'black' --define a 'black' ethnic group
  WHEN adm.ethnicity IN
    (
      'HISPANIC OR LATINO' 
    , 'HISPANIC/LATINO - PUERTO RICAN' 
    , 'HISPANIC/LATINO - DOMINICAN' 
    , 'HISPANIC/LATINO - GUATEMALAN' 
    , 'HISPANIC/LATINO - CUBAN' 
    , 'HISPANIC/LATINO - SALVADORAN' 
    , 'HISPANIC/LATINO - CENTRAL AMERICAN (OTHER)' 
    , 'HISPANIC/LATINO - MEXICAN' 
    , 'HISPANIC/LATINO - COLOMBIAN' 
    , 'HISPANIC/LATINO - HONDURAN' 
  ) THEN 'hispanic' --define an 'hispanic' ethnic group
  WHEN adm.ethnicity IN
  (
      'ASIAN'
    , 'ASIAN - CHINESE' 
    , 'ASIAN - ASIAN INDIAN'
    , 'ASIAN - VIETNAMESE' 
    , 'ASIAN - FILIPINO' 
    , 'ASIAN - CAMBODIAN' 
    , 'ASIAN - OTHER' 
    , 'ASIAN - KOREAN' 
    , 'ASIAN - JAPANESE' 
    , 'ASIAN - THAI' 
  ) THEN 'asian' --define an 'asian' ethnic group
  WHEN adm.ethnicity IN
  (
      'UNKNOWN/NOT SPECIFIED' 
    , 'UNABLE TO OBTAIN' 
    , 'PATIENT DECLINED TO ANSWER' 
  ) THEN 'unknown' --define an unknown group
  ELSE 'other' 
  END AS ethnicity_grouped
  
-- First 24 hour lab results (lfd prefix = "lab first day")
, lfd.BICARBONATE_Min
, lfd.BICARBONATE_Mean
, lfd.BICARBONATE_Max
, lfd.CREATININE_Min
, lfd.CREATININE_Mean
, lfd.CREATININE_Max
, lfd.CHLORIDE_Min
, lfd.CHLORIDE_Mean
, lfd.CHLORIDE_Max
, lfd.HEMOGLOBIN_Min
, lfd.HEMOGLOBIN_Mean
, lfd.HEMOGLOBIN_Max
, lfd.PLATELET_Min
, lfd.PLATELET_Mean
, lfd.PLATELET_Max
, lfd.POTASSIUM_Min
, lfd.POTASSIUM_Mean
, lfd.POTASSIUM_Max
, lfd.PTT_Min
, lfd.PTT_Mean
, lfd.PTT_Max
, lfd.INR_Min
, lfd.INR_Mean
, lfd.INR_Max
, lfd.PT_Min
, lfd.PT_Mean
, lfd.PT_Max
, lfd.BUN_Min
, lfd.BUN_Mean
, lfd.BUN_Max
, lfd.WBC_Min
, lfd.WBC_Mean
, lfd.WBC_Max

FROM icustays ie

INNER JOIN admissions adm
	ON ie.hadm_id = adm.hadm_id

INNER JOIN patients pat
	ON ie.subject_id = pat.subject_id

INNER JOIN(
	SELECT pvt.subject_id, pvt.hadm_id, pvt.icustay_id
	--Summarise each patient's lab tests over the first 24 hours for a given ICU admission
	, min(CASE WHEN lab_label = 'BICARBONATE' THEN valuenum ELSE null END) AS BICARBONATE_Min
	, avg(CASE WHEN lab_label = 'BICARBONATE' THEN valuenum ELSE null END) AS BICARBONATE_Mean
	, max(CASE WHEN lab_label = 'BICARBONATE' THEN valuenum ELSE null END) AS BICARBONATE_Max
	, min(CASE WHEN lab_label = 'CREATININE' THEN valuenum ELSE null END) AS CREATININE_Min
	, avg(CASE WHEN lab_label = 'CREATININE' THEN valuenum ELSE null END) AS CREATININE_Mean
	, max(CASE WHEN lab_label = 'CREATININE' THEN valuenum ELSE null END) AS CREATININE_Max
	, min(CASE WHEN lab_label = 'CHLORIDE' THEN valuenum ELSE null END) AS CHLORIDE_Min
	, avg(CASE WHEN lab_label = 'CHLORIDE' THEN valuenum ELSE null END) AS CHLORIDE_Mean
	, max(CASE WHEN lab_label = 'CHLORIDE' THEN valuenum ELSE null END) AS CHLORIDE_Max
	, min(CASE WHEN lab_label = 'HEMOGLOBIN' THEN valuenum ELSE null END) AS HEMOGLOBIN_Min
	, avg(CASE WHEN lab_label = 'HEMOGLOBIN' THEN valuenum ELSE null END) AS HEMOGLOBIN_Mean
	, max(CASE WHEN lab_label = 'HEMOGLOBIN' THEN valuenum ELSE null END) AS HEMOGLOBIN_Max
	, min(CASE WHEN lab_label = 'PLATELET' THEN valuenum ELSE null END) AS PLATELET_Min
	, avg(CASE WHEN lab_label = 'PLATELET' THEN valuenum ELSE null END) AS PLATELET_Mean
	, max(CASE WHEN lab_label = 'PLATELET' THEN valuenum ELSE null END) AS PLATELET_Max
	, min(CASE WHEN lab_label = 'POTASSIUM' THEN valuenum ELSE null END) AS POTASSIUM_Min
	, avg(CASE WHEN lab_label = 'POTASSIUM' THEN valuenum ELSE null END) AS POTASSIUM_Mean
	, max(CASE WHEN lab_label = 'POTASSIUM' THEN valuenum ELSE null END) AS POTASSIUM_Max
	, min(CASE WHEN lab_label = 'PTT' THEN valuenum ELSE null END) AS PTT_Min
	, avg(CASE WHEN lab_label = 'PTT' THEN valuenum ELSE null END) AS PTT_Mean
	, max(CASE WHEN lab_label = 'PTT' THEN valuenum ELSE null END) AS PTT_Max
	, min(CASE WHEN lab_label = 'INR' THEN valuenum ELSE null END) AS INR_Min
	, avg(CASE WHEN lab_label = 'INR' THEN valuenum ELSE null END) AS INR_Mean
	, max(CASE WHEN lab_label = 'INR' THEN valuenum ELSE null END) AS INR_Max
	, min(CASE WHEN lab_label = 'PT' THEN valuenum ELSE null END) AS PT_Min
	, avg(CASE WHEN lab_label = 'PT' THEN valuenum ELSE null END) AS PT_Mean
	, max(CASE WHEN lab_label = 'PT' THEN valuenum ELSE null END) AS PT_Max
	, min(CASE WHEN lab_label = 'BUN' THEN valuenum ELSE null end) AS BUN_Min
	, avg(CASE WHEN lab_label = 'BUN' THEN valuenum ELSE null end) AS BUN_Mean
	, max(CASE WHEN lab_label = 'BUN' THEN valuenum ELSE null end) AS BUN_Max
	, min(CASE WHEN lab_label = 'WBC' THEN valuenum ELSE null end) AS WBC_Min
	, avg(CASE WHEN lab_label = 'WBC' THEN valuenum ELSE null end) AS WBC_Mean
	, max(CASE WHEN lab_label = 'WBC' THEN valuenum ELSE null end) AS WBC_Max

	FROM( -- begin query that extracts the data:
	     SELECT ie.subject_id, ie.hadm_id, ie.icustay_id
	     -- combine multiple ITEMIDs containing the same data:
	     , CASE 
		WHEN itemid = 50882 THEN 'BICARBONATE'
		WHEN itemid = 50912 THEN 'CREATININE'
		WHEN itemid = 50806 THEN 'CHLORIDE'
		WHEN itemid = 50902 THEN 'CHLORIDE'
		WHEN itemid = 50811 THEN 'HEMOGLOBIN'
		WHEN itemid = 51222 THEN 'HEMOGLOBIN'
		WHEN itemid = 51265 THEN 'PLATELET'
		WHEN itemid = 50822 THEN 'POTASSIUM'
		WHEN itemid = 50971 THEN 'POTASSIUM'
		WHEN itemid = 51275 THEN 'PTT'
		WHEN itemid = 51237 THEN 'INR'
		WHEN itemid = 51274 THEN 'PT'
		WHEN itemid = 51006 THEN 'BUN'
		WHEN itemid = 51300 THEN 'WBC'
		WHEN itemid = 51301 THEN 'WBC'
	    ELSE null
            END AS lab_label

            -- Conduct sanity checks on the values
            , CASE 
              WHEN itemid = 50882 and valuenum > 10000 THEN null -- mEq/L 'BICARBONATE'
	      WHEN itemid = 50806 and valuenum > 10000 THEN null -- mEq/L 'CHLORIDE'
	      WHEN itemid = 50902 and valuenum > 10000 THEN null -- mEq/L 'CHLORIDE'
              WHEN itemid = 50912 and valuenum >   150 THEN null -- mg/dL 'CREATININE'
              WHEN itemid = 50811 and valuenum >    50 THEN null -- g/dL 'HEMOGLOBIN'
              WHEN itemid = 51222 and valuenum >    50 THEN null -- g/dL 'HEMOGLOBIN'
              WHEN itemid = 51265 and valuenum > 10000 THEN null -- K/uL 'PLATELET'
              WHEN itemid = 50822 and valuenum >    30 THEN null -- mEq/L 'POTASSIUM'
              WHEN itemid = 50971 and valuenum >    30 THEN null -- mEq/L 'POTASSIUM'
              WHEN itemid = 51275 and valuenum >   150 THEN null -- sec 'PTT'
              WHEN itemid = 51237 and valuenum >    50 THEN null -- 'INR'
              WHEN itemid = 51274 and valuenum >   150 THEN null -- sec 'PT'
              WHEN itemid = 51006 and valuenum >   300 THEN null -- 'BUN'
              WHEN itemid = 51300 and valuenum >  1000 THEN null -- 'WBC'
              WHEN itemid = 51301 and valuenum >  1000 THEN null -- 'WBC' 
            ELSE le.valuenum
            END AS valuenum

            FROM icustays ie
            LEFT JOIN labevents le
		ON le.subject_id = ie.subject_id 
		AND le.hadm_id = ie.hadm_id
		AND le.charttime BETWEEN (ie.intime - interval '6' hour) AND (ie.intime + interval '1' day) --extract the lab events occuring within 1 day (minus a small buffer to capture any 'pre ICU lab tests')
                AND le.ITEMID in --Only extract the lab events we are interested in for this study:
                (
                50882, -- BICARBONATE 
		50912, -- CREATININE 
		50902, -- CHLORIDE 
		50806, -- CHLORIDE, WHOLE BLOOD 
		51222, -- HEMOGLOBIN 
		50811, -- HEMOGLOBIN
		51265, -- PLATELET COUNT 
		50971, -- POTASSIUM
		50822, -- POTASSIUM, WHOLE BLOOD
		51275, -- PARTIAL THROMBOPLASTIN TIME (PTT) 
		51237, -- INTERNATIONAL NORMALIZED RATIO (INR) 
		51274, -- PROTHROMBIN TIME (PT) 
		51006, -- BLOOD UREA NITROGEN (BUN)
		51301, -- WHITE BLOOD CELLS (WBC)
		51300  -- WBC COUNT 
                )
                AND valuenum IS NOT null AND valuenum > 0 -- lab values cannot be 0 and cannot be negative
	    ) pvt
	GROUP BY pvt.subject_id, pvt.hadm_id, pvt.icustay_id --ensure summarises of lab values arefor each patient/hospitalisation/ICU event 
	) lfd --prefix with lfd (lab first day)
	ON (ie.subject_id = lfd.subject_id AND 
	    ie.hadm_id = lfd.hadm_id AND 
            ie.icustay_id = lfd.icustay_id)
-- For visual convenience of the data order the rows by subject and the admission times
ORDER BY ie.subject_id, adm.admittime, ie.intime
)
SELECT *
FROM cohort
-- Apply inclusion/exclusion criteria to the extracted cohort:
WHERE age >= 18 --only extract patients aged over 18 years (i.e. exclude neonatal)
  AND icu_los_days > 1 --exclude any ICU admission with a LOS less than 24 hours
  AND hospstay_seq = 1 --only include a patient's first admission to hospital (i.e. avoid repeated hospitalisations)
  AND icustay_seq = 1 --only consider the lab results data from a patient's first ICU admission within a given hospitalisation
;