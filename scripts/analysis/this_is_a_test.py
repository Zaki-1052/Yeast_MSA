import TC_Visualization as from datetime import datetime

from pytz import timezone


def timezone_convertion(input_datetime: datetime, input_time_zone: str):

    input_time_zone = timezone(input_time_zone)

    return input_datetime.astimezone(input_time_zone)



def get_date(date_str: str, time_str: str):

    # convert strings in format `date_str` = YYYY-MM-DD and `time_str` = HH:MM:SS into datetime

    date = list(map(int, date_str.split("-")))

    time = list(map(int, time_str.split(":")))

    return datetime(*date, *time)


TCV.say_helo()