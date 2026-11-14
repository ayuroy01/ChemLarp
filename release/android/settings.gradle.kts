// Android app: the ChemLarp web app (WebAssembly engine and UI) in a WebView. Built by the packages workflow
// (scripts/build_android.sh copies the built web app into app/src/main/assets/web first).
pluginManagement {
    repositories {
        google()
        mavenCentral()
        gradlePluginPortal()
    }
}
dependencyResolutionManagement {
    repositoriesMode.set(RepositoriesMode.FAIL_ON_PROJECT_REPOS)
    repositories {
        google()
        mavenCentral()
    }
}
rootProject.name = "ChemLarp"
include(":app")
