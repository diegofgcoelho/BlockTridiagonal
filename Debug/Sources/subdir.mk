################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Sources/DenseMatrix.cpp \
../Sources/SBlockPent.cpp \
../Sources/SBlockTrid.cpp \
../Sources/Support.cpp \
../Sources/SymMatrix.cpp \
../Sources/main.cpp 

OBJS += \
./Sources/DenseMatrix.o \
./Sources/SBlockPent.o \
./Sources/SBlockTrid.o \
./Sources/Support.o \
./Sources/SymMatrix.o \
./Sources/main.o 

CPP_DEPS += \
./Sources/DenseMatrix.d \
./Sources/SBlockPent.d \
./Sources/SBlockTrid.d \
./Sources/Support.d \
./Sources/SymMatrix.d \
./Sources/main.d 


# Each subdirectory must supply rules for building sources it contributes
Sources/%.o: ../Sources/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


