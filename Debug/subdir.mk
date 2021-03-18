################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../rinex_reader.c 

OBJS += \
./rinex_reader.o 

C_DEPS += \
./rinex_reader.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C Compiler'
	gcc -I"C:\Users\ANMT\eclipse-workspace\Garlic_IBPL_ref\garlic" -I"C:\Users\ANMT\eclipse-workspace\Garlic_IBPL_ref\ARAIM" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


