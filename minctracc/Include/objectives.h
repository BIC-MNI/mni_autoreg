float xcorr_objective(VIO_Volume d1,
                             VIO_Volume d2,
                             VIO_Volume m1,
                             VIO_Volume m2, 
                             Arg_Data *globals);

float zscore_objective(VIO_Volume d1,
                              VIO_Volume d2,
                              VIO_Volume m1,
                              VIO_Volume m2, 
                              Arg_Data *globals);

float vr_objective(VIO_Volume d1,
                          VIO_Volume d2,
                          VIO_Volume m1,
                          VIO_Volume m2, 
                          Arg_Data *globals);

float ssc_objective(VIO_Volume d1,
                           VIO_Volume d2,
                           VIO_Volume m1,
                           VIO_Volume m2, 
                           Arg_Data *globals);

float mutual_information_objective(VIO_Volume d1,
                           VIO_Volume d2,
                           VIO_Volume m1,
                           VIO_Volume m2, 
                           Arg_Data *globals);

float normalized_mutual_information_objective(VIO_Volume d1,
                           VIO_Volume d2,
                           VIO_Volume m1,
                           VIO_Volume m2, 
                           Arg_Data *globals);


