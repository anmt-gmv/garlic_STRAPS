################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../ARAIM/matrix_vect/cofactor_matrix.c \
../ARAIM/matrix_vect/copy_array_double.c \
../ARAIM/matrix_vect/copy_vect_double.c \
../ARAIM/matrix_vect/copy_vect_int.c \
../ARAIM/matrix_vect/det.c \
../ARAIM/matrix_vect/drop_matrix_col.c \
../ARAIM/matrix_vect/drop_matrix_row.c \
../ARAIM/matrix_vect/drop_vector_el.c \
../ARAIM/matrix_vect/drop_vector_el_int.c \
../ARAIM/matrix_vect/float_to_double.c \
../ARAIM/matrix_vect/init_array_double.c \
../ARAIM/matrix_vect/inv.c \
../ARAIM/matrix_vect/matrix_diff.c \
../ARAIM/matrix_vect/matrix_mul.c \
../ARAIM/matrix_vect/matrix_sum.c \
../ARAIM/matrix_vect/minor_matrix.c \
../ARAIM/matrix_vect/rank1up.c \
../ARAIM/matrix_vect/rank1up_col.c \
../ARAIM/matrix_vect/sum_elements_vect.c \
../ARAIM/matrix_vect/transpose.c 

OBJS += \
./ARAIM/matrix_vect/cofactor_matrix.o \
./ARAIM/matrix_vect/copy_array_double.o \
./ARAIM/matrix_vect/copy_vect_double.o \
./ARAIM/matrix_vect/copy_vect_int.o \
./ARAIM/matrix_vect/det.o \
./ARAIM/matrix_vect/drop_matrix_col.o \
./ARAIM/matrix_vect/drop_matrix_row.o \
./ARAIM/matrix_vect/drop_vector_el.o \
./ARAIM/matrix_vect/drop_vector_el_int.o \
./ARAIM/matrix_vect/float_to_double.o \
./ARAIM/matrix_vect/init_array_double.o \
./ARAIM/matrix_vect/inv.o \
./ARAIM/matrix_vect/matrix_diff.o \
./ARAIM/matrix_vect/matrix_mul.o \
./ARAIM/matrix_vect/matrix_sum.o \
./ARAIM/matrix_vect/minor_matrix.o \
./ARAIM/matrix_vect/rank1up.o \
./ARAIM/matrix_vect/rank1up_col.o \
./ARAIM/matrix_vect/sum_elements_vect.o \
./ARAIM/matrix_vect/transpose.o 

C_DEPS += \
./ARAIM/matrix_vect/cofactor_matrix.d \
./ARAIM/matrix_vect/copy_array_double.d \
./ARAIM/matrix_vect/copy_vect_double.d \
./ARAIM/matrix_vect/copy_vect_int.d \
./ARAIM/matrix_vect/det.d \
./ARAIM/matrix_vect/drop_matrix_col.d \
./ARAIM/matrix_vect/drop_matrix_row.d \
./ARAIM/matrix_vect/drop_vector_el.d \
./ARAIM/matrix_vect/drop_vector_el_int.d \
./ARAIM/matrix_vect/float_to_double.d \
./ARAIM/matrix_vect/init_array_double.d \
./ARAIM/matrix_vect/inv.d \
./ARAIM/matrix_vect/matrix_diff.d \
./ARAIM/matrix_vect/matrix_mul.d \
./ARAIM/matrix_vect/matrix_sum.d \
./ARAIM/matrix_vect/minor_matrix.d \
./ARAIM/matrix_vect/rank1up.d \
./ARAIM/matrix_vect/rank1up_col.d \
./ARAIM/matrix_vect/sum_elements_vect.d \
./ARAIM/matrix_vect/transpose.d 


# Each subdirectory must supply rules for building sources it contributes
ARAIM/matrix_vect/%.o: ../ARAIM/matrix_vect/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I"C:\Users\ANMT\eclipse-workspace\Garlic_IBPL_ref\garlic" -I"C:\Users\ANMT\eclipse-workspace\Garlic_IBPL_ref\ARAIM" -I"C:\MinGW\include" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


