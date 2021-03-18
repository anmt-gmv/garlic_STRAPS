################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../ARAIM/ECEF_to_ENU.c \
../ARAIM/check_new_sat_in.c \
../ARAIM/from_rinex_structures.c \
../ARAIM/load_ARAIM_parameters.c 

OBJS += \
./ARAIM/ECEF_to_ENU.o \
./ARAIM/check_new_sat_in.o \
./ARAIM/from_rinex_structures.o \
./ARAIM/load_ARAIM_parameters.o 

C_DEPS += \
./ARAIM/ECEF_to_ENU.d \
./ARAIM/check_new_sat_in.d \
./ARAIM/from_rinex_structures.d \
./ARAIM/load_ARAIM_parameters.d 


# Each subdirectory must supply rules for building sources it contributes
ARAIM/%.o: ../ARAIM/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C Compiler'
	gcc -I"C:\Users\ANMT\eclipse-workspace\Garlic_IBPL_ref\garlic" -I"C:\Users\ANMT\eclipse-workspace\Garlic_IBPL_ref\ARAIM" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


