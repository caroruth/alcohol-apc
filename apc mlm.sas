
data alc;
set apc.alcohol deaths;
where single_year_age <= 74;
cohort = year - single_year_age;
log_pop = log(population);
run;


proc glimmix data = alc;
class year cohort single_year_age (ref = '50');
model deaths = single_year_age/solution CL dist = poisson link = log offset = log_pop;
random year cohort/solution CL;
covtest GLM/WALD;
NLOPTIONS TECHNIQUE = NRRIDG;
run;
