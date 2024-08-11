% This is range of dates in file names!
first_date = "06-11-2023"
last_date = "01-14-2024"

first_date_form = datetime(first_date, "InputFormat","MM-dd-yyyy")
last_date_form = datetime(last_date, "InputFormat","MM-dd-yyyy")

for horizon = 1:4
    for date = first_date_form:7:(last_date_form-horizon*7)
    text_date = char(date, "MM-dd-yyyy")
    disp(text_date)
    plotForecast_subepidemicFramework(1, text_date, horizon, -1)
    end
end