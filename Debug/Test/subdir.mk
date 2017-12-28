################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Test/DenseMatrix_test.cpp \
../Test/SBlockPent_test.cpp \
../Test/SBlockTrid_test.cpp \
../Test/SymMatrix_test.cpp 

OBJS += \
./Test/DenseMatrix_test.o \
./Test/SBlockPent_test.o \
./Test/SBlockTrid_test.o \
./Test/SymMatrix_test.o 

CPP_DEPS += \
./Test/DenseMatrix_test.d \
./Test/SBlockPent_test.d \
./Test/SBlockTrid_test.d \
./Test/SymMatrix_test.d 


# Each subdirectory must supply rules for building sources it contributes
Test/%.o: ../Test/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


