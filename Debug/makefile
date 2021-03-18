################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include garlic/subdir.mk
-include ARAIM/matrix_vect/subdir.mk
-include ARAIM/hard_debugging/subdir.mk
-include ARAIM/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: Garlic_IBPL_ref.exe

dependents:
	-cd C:/Users/ANMT/eclipse-workspace/Garlic_IBPL_ref/Release && $(MAKE) clean all

# Tool invocations
Garlic_IBPL_ref.exe: $(OBJS) $(USER_OBJS) C:/Users/ANMT/eclipse-workspace/Garlic_IBPL_ref/Release/Garlic_IBPL_ref.exe
	@echo 'Building target: $@'
	@echo 'Invoking: Cygwin C Linker'
	gcc  -o "Garlic_IBPL_ref.exe" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(EXECUTABLES)$(OBJS)$(C_DEPS) Garlic_IBPL_ref.exe
	-@echo ' '

.PHONY: all clean dependents
C:/Users/ANMT/eclipse-workspace/Garlic_IBPL_ref/Release/Garlic_IBPL_ref.exe:

-include ../makefile.targets
