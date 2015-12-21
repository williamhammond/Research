################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Heart.cpp \
../Matrix.cpp \
../external.cpp \
../main.cpp \
../matEg.cpp 

OBJS += \
./Heart.o \
./Matrix.o \
./external.o \
./main.o \
./matEg.o 

CPP_DEPS += \
./Heart.d \
./Matrix.d \
./external.d \
./main.d \
./matEg.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/Applications/matlab_R2011b.app/extern/include/ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


