function no_nan=clearnan(have_nan)
%%ȥ��nanֵ
no_nan=have_nan(~isnan(have_nan));
end