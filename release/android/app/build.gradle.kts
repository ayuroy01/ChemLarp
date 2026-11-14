plugins {
    id("com.android.application")
}

android {
    namespace = "io.github.ayuroy01.chemlarp"
    compileSdk = 35

    defaultConfig {
        applicationId = "io.github.ayuroy01.chemlarp"
        minSdk = 26          // Android 8: WebView with WebAssembly and workers
        targetSdk = 35
        versionCode = (project.findProperty("chemlarpVersionCode") as String? ?: "1").toInt()
        versionName = project.findProperty("chemlarpVersion") as String? ?: "0.0.0"
    }

    buildTypes {
        release {
            isMinifyEnabled = false
        }
    }

    compileOptions {
        sourceCompatibility = JavaVersion.VERSION_17
        targetCompatibility = JavaVersion.VERSION_17
    }

    androidResources {
        noCompress += listOf("wasm")  // read straight from the APK
    }
}
