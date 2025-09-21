Herramienta bioinformática que automatiza el análisis de experimentos transcriptómicos de ArrayExpress: descarga, preprocesa y transforma los datos de expresión génica para hacerlos compatibles con el algoritmo de triclustering TriGen. Incluye visualizaciones de resultados e integración con PantherDB para análisis de enriquecimiento funcional (Gene Ontology). 
GitHub.

Tecnologías: R (app tipo Shiny) con lógica auxiliar en helpers.R. El repositorio contiene además material de soporte (presentación, memoria y vídeo).

Funcionalidades
- Ingesta desde ArrayExpress: descarga de datos de expresión a partir de un identificador/estudio.
- Preprocesado automático: limpieza, selección/transformación de genes y normalización.
- Preparación para TriGen: generación del formato de entrada esperado por el algoritmo.
- Triclustering: ejecución del flujo analítico para detectar triclusters (genes × condiciones × tiempos).
- Análisis funcional (PantherDB/GO): enriquecimiento de términos GO para interpretar los triclusters.
- Visualización: gráficos para explorar calidad, patrones y resultados del triclustering.
(Las capacidades se resumen desde la descripción del proyecto en el repositorio.)

Flujo de trabajo
1. Seleccionar estudio de ArrayExpress (ID o metadatos).
2. Descargar y preprocesar la matriz de expresión.
3. Formatear para TriGen (estructuras/archivos intermedios).
4. Ejecutar triclustering y explorar visualizaciones.
5. Lanzar enriquecimiento GO en PantherDB y revisar informes.

Este repositorio incluye memoria y presentación del TFG, así como un vídeo de la herramienta en funcionamiento, que contextualizan objetivos, decisiones técnicas y resultados obtenidos.
