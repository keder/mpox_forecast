first_date = "06-11-2023"
last_date = "01-14-2024"

first_date_form = datetime(first_date, "InputFormat","MM-dd-yyyy")
last_date_form = datetime(last_date, "InputFormat","MM-dd-yyyy")

pool = parpool("Processes")
% pool = parpool(36)

i=1;
for date = first_date_form:7:last_date_form
    text_date = char(date, "MM-dd-yyyy")
    disp(text_date)
    % Run_Fit_subepidemicFramework(1, text_date)
    f(i) = parfeval(@Run_Fit_subepidemicFramework, 0, 1, text_date)
    i = i + 1;
end
wait(f)
delete(pool)
