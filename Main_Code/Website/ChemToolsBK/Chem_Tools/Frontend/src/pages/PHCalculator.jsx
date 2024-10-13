import React, { useState } from 'react';
import { Link } from 'react-router-dom';
import '../App.css';

const PHCalculator = () => {
  const [concentracao, setConcentracao] = useState('');
  const [pH, setPH] = useState('');

  const handleSubmit = (e) => {
    e.preventDefault();
    // Funcionalidade para calcular o pH será implementada aqui
    setPH('O pH calculado aparecerá aqui');
  };

  return (
    <div className="ferramentas-quimica">
      <h1>Calculadora de pH</h1>
      <p className="description">Ferramenta essencial para cálculos de pH em soluções</p>
      
      <div className="conteudo">
        <div className="button-grid">
          <div className="category-column">
            <h2>Calculadora</h2>
            <form onSubmit={handleSubmit}>
              <div>
                <label htmlFor="concentracao">Digite a concentração de íons H+ (mol/L):</label>
                <input
                  type="number"
                  id="concentracao"
                  value={concentracao}
                  onChange={(e) => setConcentracao(e.target.value)}
                  placeholder="ex., 1,0 x 10^-7"
                  className="aba input-equation"
                />
              </div>
              <button type="submit" className="aba">Calcular pH</button>
            </form>
            {pH && (
              <div className="result">
                <h3>pH:</h3>
                <p>{pH}</p>
              </div>
            )}
          </div>
          <div className="category-column">
            <h2>Descrição do Problema</h2>
            <p>
              O pH é uma medida da acidez ou basicidade de uma solução aquosa. Esta ferramenta
              ajuda você a calcular o pH de uma solução com base na concentração de íons H+,
              fornecendo informações cruciais para experimentos e análises químicas.
            </p>
            <p>
              Para usar, insira a concentração de íons H+ em mol/L no campo de entrada e
              clique em "Calcular pH". A ferramenta irá determinar o pH da solução, que varia
              de 0 (muito ácido) a 14 (muito básico), com 7 sendo neutro.
            </p>
          </div>
        </div>
      </div>
      <Link to="/" className="aba">Voltar para a Página Inicial</Link>
    </div>
  );
};

export default PHCalculator;
