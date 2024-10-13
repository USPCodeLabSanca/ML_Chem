import React, { useState } from 'react';
import { Link } from 'react-router-dom';
import '../App.css';

const CurvaTitulacao = () => {
  const [acido, setAcido] = useState('');
  const [base, setBase] = useState('');
  const [concentracaoAcido, setConcentracaoAcido] = useState('');
  const [concentracaoBase, setConcentracaoBase] = useState('');
  const [volumeAcido, setVolumeAcido] = useState('');
  const [dadosCurva, setDadosCurva] = useState(null);

  const handleSubmit = (e) => {
    e.preventDefault();
    // Funcionalidade para calcular a curva de titulação será implementada aqui
    setDadosCurva('Os dados da curva de titulação serão exibidos aqui');
  };

  return (
    <div className="ferramentas-quimica">
      <h1>Calculadora de Curva de Titulação</h1>
      <p className="description">Ferramenta essencial para calcular curvas de titulação</p>
      
      <div className="conteudo">
        <div className="button-grid">
          <div className="category-column">
            <h2>Calculadora</h2>
            <form onSubmit={handleSubmit}>
              <div>
                <label htmlFor="acido">Ácido:</label>
                <input
                  type="text"
                  id="acido"
                  value={acido}
                  onChange={(e) => setAcido(e.target.value)}
                  placeholder="ex., HCl"
                  className="aba input-equation"
                />
              </div>
              <div>
                <label htmlFor="base">Base:</label>
                <input
                  type="text"
                  id="base"
                  value={base}
                  onChange={(e) => setBase(e.target.value)}
                  placeholder="ex., NaOH"
                  className="aba input-equation"
                />
              </div>
              <div>
                <label htmlFor="concentracaoAcido">Concentração do Ácido (M):</label>
                <input
                  type="number"
                  id="concentracaoAcido"
                  value={concentracaoAcido}
                  onChange={(e) => setConcentracaoAcido(e.target.value)}
                  placeholder="ex., 0,1"
                  className="aba input-equation"
                />
              </div>
              <div>
                <label htmlFor="concentracaoBase">Concentração da Base (M):</label>
                <input
                  type="number"
                  id="concentracaoBase"
                  value={concentracaoBase}
                  onChange={(e) => setConcentracaoBase(e.target.value)}
                  placeholder="ex., 0,1"
                  className="aba input-equation"
                />
              </div>
              <div>
                <label htmlFor="volumeAcido">Volume Inicial do Ácido (mL):</label>
                <input
                  type="number"
                  id="volumeAcido"
                  value={volumeAcido}
                  onChange={(e) => setVolumeAcido(e.target.value)}
                  placeholder="ex., 50"
                  className="aba input-equation"
                />
              </div>
              <button type="submit" className="aba">Calcular Curva de Titulação</button>
            </form>
            {dadosCurva && (
              <div className="result">
                <h3>Curva de Titulação:</h3>
                <p>{dadosCurva}</p>
                {/* Um componente de gráfico será adicionado aqui para exibir a curva de titulação */}
              </div>
            )}
          </div>
          <div className="category-column">
            <h2>Descrição do Problema</h2>
            <p>
              A curva de titulação é uma representação gráfica da mudança de pH durante uma titulação
              ácido-base. Esta ferramenta ajuda você a calcular e visualizar a curva de titulação para
              uma dada combinação de ácido e base, permitindo uma melhor compreensão do processo de
              neutralização e dos pontos críticos como o ponto de equivalência.
            </p>
            <p>
              Para usar, insira as informações sobre o ácido e a base que você está utilizando,
              incluindo suas concentrações e o volume inicial do ácido. Clique em "Calcular Curva de Titulação"
              e a ferramenta irá gerar os dados necessários para plotar a curva de titulação.
            </p>
          </div>
        </div>
      </div>
      <Link to="/" className="aba">Voltar para a Página Inicial</Link>
    </div>
  );
};

export default CurvaTitulacao;
