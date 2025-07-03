#define inputPin 13
#define inputPin1 0
#define inputPin2 0
#define outputPin0 15
#define outputPin 16
#define outputPin1 0
#define outputPin2 0
#define outputPin3 0
#define timeOUT 0
#define firePin 4
#define STAGES 2

volatile bool fire = false;
volatile bool triggerDetected = false;
volatile uint32_t timeStamp = 0;

bool measuring = false;
uint32_t startTime = 0;
uint32_t elapsedTime = 0;
float velocity = 0.0;

// Settings
const float projectileLength_m = 0.06; 
const float distanceToTravel_m = 0.02;
volatile int count = 0;

void IRAM_ATTR fireInterrupt() {
  fire = true;
}


void IRAM_ATTR handleInterrupt() {
  if (!measuring) {
    startTime = micros();
    measuring = true;
  }
}

void IRAM_ATTR handleFallInterrupt() {
    timeStamp = micros();
    triggerDetected = true;
}


void setup() {
  Serial.begin(115200);
  
  pinMode(firePin, INPUT_PULLDOWN);
  pinMode(inputPin, INPUT_PULLUP);
  //pinMode(inputPin1, INPUT_PULLUP);
  pinMode(outputPin, OUTPUT);
  pinMode(outputPin0, OUTPUT);

  digitalWrite(outputPin, LOW);
  digitalWrite(outputPin0, LOW);
  
  attachInterrupt(digitalPinToInterrupt(firePin), fireInterrupt, RISING);
  attachInterrupt(digitalPinToInterrupt(inputPin), handleInterrupt, RISING);
   attachInterrupt(digitalPinToInterrupt(inputPin), handleFallInterrupt, FALLING);
  //attachInterrupt(digitalPinToInterrupt(inputPin1), handleInterrupt, CHANGE);
}



void loop() {
  if (fire) {
    stage1();
  }
  if (triggerDetected) {
    //stagex();
  }
}


void stage1() {
  digitalWrite(outputPin0, HIGH);
  digitalWrite(timeOUT, HIGH);
  delay(1);
  digitalWrite(timeOUT, LOW);
  digitalWrite(outputPin0, LOW);
  Serial.println("firing");

  fire = false;
}

void stagex() {
  elapsedTime = timeStamp - startTime; // microseconds
  float elapsedTime_s = elapsedTime / 1e6; // convert to seconds
  velocity = projectileLength_m / elapsedTime_s; // meters per second
  float waitTime_s = distanceToTravel_m / velocity; // Now wait for projectile to travel distanceToTravel_m
  uint32_t waitTime_ms = (uint32_t)(waitTime_s * 1000000.0); // convert to micros
  delayMicroseconds(waitTime_ms); // blocking wait
  digitalWrite(outputPin, HIGH);
  digitalWrite(timeOUT, HIGH);
  delay(1);
  digitalWrite(timeOUT, LOW);
  digitalWrite(outputPin, LOW);
  count = count + 1;
  //Serial.println("fire stage");
  measuring = false; // Reset for next measurement
  triggerDetected = false;
}
