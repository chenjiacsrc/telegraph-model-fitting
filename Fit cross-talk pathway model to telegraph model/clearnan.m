function no_nan=clearnan(have_nan)
%%È¥µônanÖµ
no_nan=have_nan(~isnan(have_nan));
end