
# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/fista.cpp \
../src/main.cpp \
../src/recon.cpp 

OBJS += \
./src/fista.o \
./src/main.o \
./src/recon.o 

CPP_DEPS += \
./src/fista.d \
./src/main.d \
./src/recon.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -Ifftw3f -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


