/*
  ESP32-S3 Laser Distance & Velocity Logger
  Reads 0-5V output from IL-1000 sensor box with IL-2000 laser sensor on GPIO7
  Includes voltage divider scaling, millisecond time display, and improved data clearing
*/

#include <WiFi.h>
#include <WebServer.h>
//#include <SPIFFS.h>

// ============ CONFIGURATION SECTION ============
//const char* ssid = "NoahFL";
//const char* password = "noah1234";

//const char* ssid = "Ronan";
//const char* password = "nanor17!";


// Sensor settings (from IL-2000 datasheet)
const float SENSOR_MIN_DISTANCE = 1000.0;   // mm (minimum measurement range)
const float SENSOR_MAX_DISTANCE = 3000.0;   // mm (maximum measurement range)
const float SENSOR_VOLTAGE_MIN = 0.0;       // V (corresponds to 1000mm)
const float SENSOR_VOLTAGE_MAX = 5.0;       // V (corresponds to 3000mm)

// Laser Voltage measurement settings
const int LASER_ADC_PIN = 8;                // GPIO pin for laser sensor input
const int ADC_RESOLUTION = 12;              // 12-bit ADC resolution
const float ADC_REF_VOLTAGE = 3.3;          // Reference voltage for ADC
const float VOLTAGE_DIVIDER_RATIO = 1.5667;  // Voltage divider ratio
const int FIRE_PIN = 14;  // GPIO 14 for digital trigger

// For voltage measurment
const float VOLTAGE_DIVIDER = 200.0;
const int ADC_PIN_1 = 6;  // First ADC pin
const int ADC_PIN_2 = 7;  // Second ADC pin

// Sampling settings
const int SAMPLE_INTERVAL_MS = 10;          // 10ms sample interval (100Hz)
const int MAX_SAMPLES = 1000;                // Maximum number of samples to store
const int VELOCITY_WINDOW = 5;              // Number of samples to use for velocity calculation

// ============ GLOBAL VARIABLES ============
WebServer server(80);

struct Laser_Sample {
  unsigned long timestamp;  // Time in milliseconds (relative to first sample)
  float voltage;            // Raw voltage reading (after divider)
  float distance;           // Calculated distance in mm
  float velocity;           // Calculated velocity in mm/s
};

Laser_Sample samples[MAX_SAMPLES];
int sampleCount = 0;
bool samplingActive = false;
bool VsamplingActive = false;
unsigned long firstSampleTime = 0;  // To track relative time

struct DAQ_Sample {
  unsigned long timestamp; 
  float voltage1;
  float voltage2;          
  //float current;           
};

DAQ_Sample daq_samples[MAX_SAMPLES];
//int daq_sampleCount = 0;
int daq_head = 0;
int daq_size = 0;

void laserTask(void *arg);
void webServerTask(void *arg);
void voltageMonitorTask(void *arg);
float mapFloat();
void handleData();
void handleLaserDownload();
//void handleClear();
void handleRoot();



// ============ SETUP FUNCTION ============
void setup() {
  Serial.begin(115200);
  //esp_log_level_set("*", ESP_LOG_VERBOSE);
  delay(50);
  
  struct WiFiCredential {
    const char* ssid;
    const char* password;
  };

    WiFiCredential networks[] = {
    {"Ronan", "nanor17!"},
    {"NoahFL", "noah1234"},
    {"saras iphone", "m33pm00p"} // Replace with actual third network
  };

  const int numNetworks = sizeof(networks) / sizeof(networks[0]);
  bool connected = false;


  //WiFi.begin(ssid, password);
  Serial.print("Connecting to WiFi");
  for (int i = 0; i < numNetworks; ++i) {
    Serial.printf("Trying to connect to SSID: %s\n", networks[i].ssid);
    WiFi.begin(networks[i].ssid, networks[i].password);

    unsigned long startAttemptTime = millis();

    // Wait up to 8 seconds for connection
    while (WiFi.status() != WL_CONNECTED && millis() - startAttemptTime < 8000) {
      delay(500);
      Serial.print(".");
    }

    if (WiFi.status() == WL_CONNECTED) {
      Serial.println("\nConnected!");
      connected = true;
      break;
    } else {
      Serial.println("\nFailed");
    }
  }

  server.on("/", handleRoot);
  server.on("/start", []() {
    VsamplingActive = true;
    server.send(200, "text/plain", "Sampling started");
  });
  server.on("/data", handleData);
  server.on("/download", handleLaserDownload);
  server.on("/fire", []() {
    digitalWrite(FIRE_PIN, HIGH);
    //delay(1000);
    delayMicroseconds(100);
    digitalWrite(FIRE_PIN, LOW);
    Serial.println("FIRE pin toggled");

    // Reset laser data buffer
    sampleCount = 0;
    firstSampleTime = millis();

    // Insert initial "fire" marker sample at index 0
    samples[0] = {0, 0.0, 0.0, 0.0};  // You could use a sentinel value if needed
    sampleCount = 1;                 // Reserve index 0

    // Start laser sampling
    samplingActive = true;

    server.send(200, "text/plain", "FIRE triggered");
  });

  server.begin();
  Serial.println("HTTP server started");

  // Split tasks between cores
  xTaskCreatePinnedToCore(laserTask, "Laser", 4096, NULL, 1, NULL, 0);
  xTaskCreatePinnedToCore(webServerTask, "WebServer", 8192, NULL, 1, NULL, 1);
  xTaskCreatePinnedToCore(voltageMonitorTask, "Voltage", 4096, NULL, 1, NULL, 1);
  
  pinMode(FIRE_PIN, OUTPUT);
  digitalWrite(FIRE_PIN, LOW);

  // Initialize SPIFFS
  //if (!SPIFFS.begin(true)) {
  //  Serial.println("SPIFFS Mount Failed");
  //}
  // Configure ADC
  analogReadResolution(ADC_RESOLUTION);

}

// ============ MAIN LOOP ============
void loop() {

}

void voltageMonitorTask(void* pvParameters) {
  analogSetPinAttenuation(ADC_PIN_1, ADC_11db);
  analogSetPinAttenuation(ADC_PIN_2, ADC_11db);

  while (true) {
    vTaskDelay(10 / portTICK_PERIOD_MS);
    if (!VsamplingActive) {
      vTaskDelay(50 / portTICK_PERIOD_MS);
      continue;
    }

    int adc1 = analogRead(ADC_PIN_1);
    int adc2 = analogRead(ADC_PIN_2);

    float v1 = (adc1 * ADC_REF_VOLTAGE) / (1 << ADC_RESOLUTION) * VOLTAGE_DIVIDER;
    float v2 = (adc2 * ADC_REF_VOLTAGE) / (1 << ADC_RESOLUTION) * VOLTAGE_DIVIDER;

    daq_samples[daq_head] = {millis(), v1, v2};
    daq_head = (daq_head + 1) % MAX_SAMPLES;
    if (daq_size < MAX_SAMPLES) daq_size++;

    vTaskDelay(100 / portTICK_PERIOD_MS);
  }
}

void webServerTask(void* pvParameters) {
  while (true) {
    //Serial.printf("Free heap: %lu\n", esp_get_free_heap_size());
    server.handleClient();
    //Serial.printf("WiFi status: %d\n", WiFi.status());
    vTaskDelay(200 / portTICK_PERIOD_MS);
  }
}

void laserTask(void* pvParameters) {
  analogSetPinAttenuation(LASER_ADC_PIN, ADC_11db);
  while (true) {
    if (!samplingActive) {
      vTaskDelay(10 / portTICK_PERIOD_MS);
      continue;
    }
    if (sampleCount >= MAX_SAMPLES) {
      samplingActive = false;
      Serial.println("Sample buffer full");
      continue;
    }
    // Read and convert voltage from laser sensor pin (with divider scaling)
    int adcValue = analogRead(LASER_ADC_PIN);
    float voltage = (adcValue * ADC_REF_VOLTAGE) / (1 << ADC_RESOLUTION) * VOLTAGE_DIVIDER_RATIO;
    
    // Calculate distance (linear conversion from voltage to mm)
    float distance = mapFloat(voltage, 
                            SENSOR_VOLTAGE_MIN * VOLTAGE_DIVIDER_RATIO, 
                            SENSOR_VOLTAGE_MAX * VOLTAGE_DIVIDER_RATIO,
                            SENSOR_MIN_DISTANCE, SENSOR_MAX_DISTANCE);
    
    // Calculate velocity (if we have previous samples)
    float velocity = 0.0;
    if (sampleCount > VELOCITY_WINDOW) {
      int oldestIndex = max(0, sampleCount - VELOCITY_WINDOW);
      float distanceChange = distance - samples[oldestIndex].distance;
      float timeChange = (samples[sampleCount - 1].timestamp - samples[oldestIndex].timestamp) / 1000.0; // in seconds
      velocity = distanceChange / timeChange;
    }
    
    // Store the sample with relative timestamp
    samples[sampleCount] = {millis() - firstSampleTime, voltage, distance, velocity};
    sampleCount++;
    
    vTaskDelay(2 / portTICK_PERIOD_MS);
  }
}

  // Helper function for floating-point mapping
float mapFloat(float x, float in_min, float in_max, float out_min, float out_max) {
  return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}

void handleData() {
  const int N = 100;
  int count = min(N, daq_size);
  int start = (daq_head - count + MAX_SAMPLES) % MAX_SAMPLES;

  String json = "{\"t\":[";
  for (int i = 0; i < count; i++) {
    int idx = (start + i) % MAX_SAMPLES;
    if (i > 0) json += ",";
    json += daq_samples[idx].timestamp;
  }

  json += "],\"v1\":[";
  for (int i = 0; i < count; i++) {
    int idx = (start + i) % MAX_SAMPLES;
    if (i > 0) json += ",";
    json += daq_samples[idx].voltage1;
  }

  json += "],\"v2\":[";
  for (int i = 0; i < count; i++) {
    int idx = (start + i) % MAX_SAMPLES;
    if (i > 0) json += ",";
    json += daq_samples[idx].voltage2;
  }

  json += "]}";
  server.send(200, "application/json", json);
}


void handleLaserDownload() {
  //String csv = "Time(ms),";
  //String filename;
  
  String csv = "Time(ms),Voltage(V),Distance(mm),Velocity(m/s)\n";
  for (int i = 0; i < sampleCount; i++) {
    csv += String(samples[i].timestamp) + "," + String(samples[i].voltage, 4) + "," + String(samples[i].distance, 2) + "," + String(samples[i].velocity, 2)+"\n";
  }
  server.send(200, "text/csv", csv);
}

/*void handleClear() {
  // Reset laser data
  sampleCount = 0;
  firstSampleTime = millis();
  
  // Reset voltage data
  daq_head = 0;
  daq_size = 0;
  
  // Send response that will trigger chart clearing
  String response = R"rawliteral({
    "v1": [],
    "v2": [],
    "t": []
  })rawliteral";
  
  server.send(200, "application/json", response);
}*/

void handleRoot() {
  String html = R"rawliteral(
  <!DOCTYPE html>
  <html>
  <head>
    <title>Voltage Logger</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
      body { font-family: Arial, sans-serif; margin: 20px; }
      button { padding: 8px 16px; margin: 5px; }
      #chart { background: #f8f8f8; }
    </style>
  </head>
  <body>
    <h1>Control Panel</h1>

    <div>
      <button onclick="startSampling()">Start</button>
      <button onclick="stopSampling()">Stop</button>
      <button onclick="downloadLaserData()">Download Laser</button> 
      <button onclick="fire()">FIRE</button>
    </div>

    <div>Status: <span id="status">Ready</span></div>
    <canvas id="chart"></canvas>
    <script>
      const ctx = document.getElementById('chart').getContext('2d');
      const chart = new Chart(ctx, {
        type: 'line',
        data: {
          labels: [],
          datasets: [
            { label: 'Voltage 1 (V)', borderColor: 'rgb(75, 192, 192)', data: [] },
            { label: 'Voltage 2 (V)', borderColor: 'rgb(255, 99, 132)', data: [] }
          ]
        },
        options: {
          animation: false,
          scales: {
            x: { title: { display: true, text: 'Time (s)' }},
            y: { title: { display: true, text: 'Voltage (V)' }}
          }
        }
      });

      function updateChart() {
        fetch('/data')
          .then(r => r.json())
          .then(data => {
            const time = data.t.map(t => (t/1000).toFixed(1));
            chart.data.labels = time;
            chart.data.datasets[0].data = data.v1;
            chart.data.datasets[1].data = data.v2;
            chart.update();
            document.getElementById('status').innerText = `Samples: ${data.v1.length}`;
          });
      }

      function startSampling() { 
        fetch('/start').then(updateChart); 
      }
      
      function stopSampling() { 
        fetch('/stop').then(updateChart); 
      }
      
      function downloadLaserData() { 
        window.location.href = '/download'; 
      }
      
      function fire() {
        fetch('/fire').then(updateChart);
      }
      
      // Update chart every 200ms
      setInterval(updateChart, 500);
      updateChart();
    </script>
  </body>
  </html>
  )rawliteral";

  server.send(200, "text/html", html);
}