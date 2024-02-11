Используя метод грамматической эволюции, в работах [1] и [2] были построены модели
времени цветения нута. В работе [3] был предложен другой подход к построению функции
регрессии - декартово генетическое программирование. В работе [4] был модифицирован
алгоритм декартова генетического программирования для учета в нем констант. В работе
[4] представлен алгоритм подбора гиперпараметров ABCDE, являющийся формой
байесовского метода оптимизации. Сами байесовские методы оптимизации дают не
точечную оценку управляющих параметров модели, а оценку их апостериорного
распределения. Из минусов — большая вычислительная сложность. Поэтому важно
вовремя остановить вычисления, с какой-то вероятностью достигнув заданной точности.
Цель этой работы – использовать алгоритм ABCDE для подбора гиперпараметров в
декартовом генетическом программировании для создания модели времени цветения,
предложив и реализовав при этом методы определения сходимости метода ABCDE. Кроме
этого, необходимо исследовать полученные решения.


1. Ageev, Andrey, Abdulkadir Aydogan, Eric Bishop-von Wettberg, Sergey V. Nuzhdin,
Maria Samsonova, and Konstantin Kozlov. “Simulation Model for Time to Flowering
with Climatic and Genetic Inputs for Wild Chickpea.” Agronomy 11, no. 7 (July 9, 2021):
1389. https://doi.org/10.3390/agronomy11071389.[Ageev_et_al-2021-
Simulation_Model_for_Time_to_Flowering_with_Climatic_and_Genetic_Inputs_for.pdf
](uploads/9910bd0f0cfc117d5c7fee2f79cbb2e7/Ageev_et_al-2021-
Simulation_Model_for_Time_to_Flowering_with_Climatic_and_Genetic_Inputs_for.pdf
)
2. Kozlov, Konstantin, Anupam Singh, Jens Berger, Eric Bishop-von Wettberg, Abdullah
Kahraman, Abdulkadir Aydogan, Douglas Cook, Sergey Nuzhdin, and Maria Samsonova.
“Non-Linear Regression Models for Time to Flowering in Wild Chickpea Combine
Genetic and Climatic Factors.” BMC Plant Biology 19, no. S2 (2019): 94.
3. Miller, Julian Francis. “Cartesian Genetic Programming: Its Status and Future.” Genetic
Programming and Evolvable Machines 21, no. 1–2 (June 2020): 129–68.
4. Михальчук В. А. Разработка модели учета генетических и географических факторов
для прогнозирования фенотипических признаков растений с помощью декартова
генетического программирования: выпускная квалификационная работа бакалавра.
5. Turner, B. M. & Sederberg, P. B. Approximate Bayesian computation with differential
evolution. Journal of Mathematical Psychology 56, 375–385 (2012).
