################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../kmeans.c \
../pqnn32c.c 

O_SRCS += \
../pqnn32.o 

OBJS += \
./kmeans.o \
./pqnn32c.o 

C_DEPS += \
./kmeans.d \
./pqnn32c.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

pqnn32c.o: ../pqnn32c.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc  -m32 -msse ../pqnn32.o ../pqnn32c.c -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"pqnn32c.d" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


