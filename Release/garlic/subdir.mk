################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../garlic/GICSRxAux.c \
../garlic/GICSRxIBPL.c \
../garlic/GICSRxInstallation.c \
../garlic/GICSRxKalman.c \
../garlic/GICSRxKalman_gnss.c \
../garlic/GICSRxMechanization.c \
../garlic/GICSRxObs.c \
../garlic/GICSRxPosition.c \
../garlic/GICSRxRINEX.c \
../garlic/GICSRxWNavSol.c \
../garlic/Indicators_Recorder.c \
../garlic/NeQuick.c \
../garlic/Residuals_Recorder.c \
../garlic/Summary_Recorder.c \
../garlic/algebra.c \
../garlic/atmsphcorr.c \
../garlic/calendar.c \
../garlic/garlic_functions.c \
../garlic/matrix.c 

OBJS += \
./garlic/GICSRxAux.o \
./garlic/GICSRxIBPL.o \
./garlic/GICSRxInstallation.o \
./garlic/GICSRxKalman.o \
./garlic/GICSRxKalman_gnss.o \
./garlic/GICSRxMechanization.o \
./garlic/GICSRxObs.o \
./garlic/GICSRxPosition.o \
./garlic/GICSRxRINEX.o \
./garlic/GICSRxWNavSol.o \
./garlic/Indicators_Recorder.o \
./garlic/NeQuick.o \
./garlic/Residuals_Recorder.o \
./garlic/Summary_Recorder.o \
./garlic/algebra.o \
./garlic/atmsphcorr.o \
./garlic/calendar.o \
./garlic/garlic_functions.o \
./garlic/matrix.o 

C_DEPS += \
./garlic/GICSRxAux.d \
./garlic/GICSRxIBPL.d \
./garlic/GICSRxInstallation.d \
./garlic/GICSRxKalman.d \
./garlic/GICSRxKalman_gnss.d \
./garlic/GICSRxMechanization.d \
./garlic/GICSRxObs.d \
./garlic/GICSRxPosition.d \
./garlic/GICSRxRINEX.d \
./garlic/GICSRxWNavSol.d \
./garlic/Indicators_Recorder.d \
./garlic/NeQuick.d \
./garlic/Residuals_Recorder.d \
./garlic/Summary_Recorder.d \
./garlic/algebra.d \
./garlic/atmsphcorr.d \
./garlic/calendar.d \
./garlic/garlic_functions.d \
./garlic/matrix.d 


# Each subdirectory must supply rules for building sources it contributes
garlic/%.o: ../garlic/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I"C:\Users\ANMT\eclipse-workspace\Garlic_IBPL_ref\garlic" -I"C:\Users\ANMT\eclipse-workspace\Garlic_IBPL_ref\ARAIM" -I"C:\MinGW\include" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


